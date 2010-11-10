#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "Common.h"
#include "EMClusteringIO.h"
#include "AffineRegistration.h"
#include "myMaths.h"
#include "EMCore.h"
#include "MeshOperations.h"
#include "Quantification.h"

#include "EMClusteringCLP.h"

int main(int argc, char* argv[])
{
	PARSE_ARGS;

	MeshType::Pointer    Trajectories;
	CenterType           Centers;
	ImageType::Pointer   subSpace;

	bool debug = 0;
	CopyFieldType copyField = {0,0,1,1,0,1};
	unsigned int ncells; //number of cells in each mesh

	// Get the input trajectories
	std::vector<std::string> allfilenames;
	if (population)
	{
		copyField.CaseName = 1;
		itksys::Glob allfiles;
		std::string globPath = TractsDir+"/*.vt*"; //vtk or vtp
		allfiles.FindFiles(globPath);
		allfilenames = allfiles.GetFiles();
		Trajectories = ReadVTKfiles(allfilenames);
	}
	else if (!trajectoriesFilename.empty())
	{
		copyField.CaseName = 0;
		Trajectories = ReadVTKfile(trajectoriesFilename.c_str());
	}
	else
	{
		std::cerr << "No input was given as a collection of trajectories to be clustered" << std::endl;
		return -1;
	}

	if (!centersFilename.empty())
    {
	   Centers = readCenterFiles(centersFilename, ncells);  //ToDO: what is ncells used for?
	}
	else
	{
		std::cerr << "No initial center was given" << std::endl;
		return -1;
	}

	if (use_atlas)
	{
		//Read FA map of the atlas
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		std::cout << "Reading " << atlasFAFilename.c_str() << " ..." <<std::endl;
		reader->SetFileName(atlasFAFilename.c_str());
		reader->Update();
		ImageType::Pointer atlasFAVolume = reader->GetOutput();
		//Read FA map of the case
		ReaderType::Pointer reader2 = ReaderType::New();
		reader2->SetFileName(subjectFAFilename.c_str());
		std::cout << "Reading " << subjectFAFilename.c_str() << " ..." <<std::endl;
		reader2->Update();
		ImageType::Pointer caseFAVolume = reader2->GetOutput();

		if (Centers.size()>0)
		{
			//Do affine registration
			TransformType::Pointer transform;
			std::cout << "Registering the atlas' scalar volume to subject's ..." << std::endl;
			transform = doSlicerFastAffineRegistration(caseFAVolume, atlasFAVolume, OutputDirectory);
			//Select and transfer the centers
			Centers = applyTransform(Centers, transform);
			CopyFieldType copyField = {0,0,0,0};
			if (!transformedCentersFilename.empty())
			{
				WriteVTKfile(Centers, transformedCentersFilename,copyField);
			}
		}
	}

	//////////////////////////////////////////////////////
	// Continue if you have valid trajectories and centers:
	//////////////////////////////////////////////////////
	unsigned int numberOfTrajectories = Trajectories->GetNumberOfCells();
	unsigned int numberOfCenters = Centers.size();

	if (numberOfTrajectories>0 && numberOfCenters>0)
	{
		AddOrientation(Trajectories);

		VariableType tractographyStepSize  = getSampleSpacing(Trajectories);
		std::cout << "Calculated tractography step size is " << tractographyStepSize << std::endl;

		// set the space to the limits of input trajectories

		VariableType spaceResolution = 3; //mm  initial value for the first iteration -- space resolution >= tractographyStepSize
	    std::cout << "Setting the initial space resolution to " << spaceResolution << std::endl;
	    samplesDistance = 3.5; //mm  initial value for the first iteration
	    // Resample initial centers:
	    bool justResample =1;
	    Centers = SmoothMeshes(Centers, samplesDistance, justResample);

		VariableType MinPost = (VariableType) 1/numberOfCenters;
		// If the number of clusters is 1:
		if (MinPost>0.5) { MinPost = 0.5;}

		VariableType MinLike = 0.1*MinLikelihoodThr;

		//EM Initialization

		ArrayType alpha, beta, MyMinLikelihoodThr,dd;
		Array2DType DissimilarityMatrix, Likelihood, Prior, Posterior; //NxK

		alpha.SetSize(numberOfCenters); alpha.fill(1);
		beta.SetSize(numberOfCenters); beta.fill(10);
		Prior.SetSize(numberOfTrajectories,numberOfCenters);
		bool havePrior = 0;
		VariableType initp = ((VariableType) 1.0) /numberOfCenters;
		Prior.fill(initp);
		VariableType cosAngleThreshold = 0.5;

		bool considerOrientation =1;
		for (int i=0; i<maxNumberOfIterations; ++i)
		{
		  std::cout<< "Iteration  " << i+1 << std::endl;
		  std::stringstream currentIteration;
		  currentIteration << i+1;

		  subSpace = getSubSpace(Trajectories, spaceResolution);

		  AddOrientation(Centers);

		  cosAngleThreshold+= i*0.1;  //make it a tighter constraint every iteration
		  cosAngleThreshold = std::min(cosAngleThreshold, (VariableType) 0.8);

		  DissimilarityMatrix = ComputeDissimilarity(Trajectories, Centers, subSpace, MaxDist, cosAngleThreshold, considerOrientation);

		  //EM Block Begins
		  Likelihood = ComputeLikelihood(DissimilarityMatrix, alpha, beta);
		  MyMinLikelihoodThr = AdjustThreshold(MinLike, alpha, beta)/(i+1);

		  if (debug)
		  {
			  std::cout<< DissimilarityMatrix << std::endl;
			  std::cout << DissimilarityMatrix.max_value() << " " << DissimilarityMatrix.min_value() << std::endl;
			  std::cout<< Likelihood << std::endl;
			  std::cout<< "MyMinLikelihoodThr = " << MyMinLikelihoodThr << std::endl;
		  }

		  MeshType::Pointer RefinedTrajectories = RefineData(Trajectories,DissimilarityMatrix,Likelihood,Prior,MyMinLikelihoodThr,havePrior);
		  std::cout << RefinedTrajectories->GetNumberOfCells() << " out of " << Trajectories->GetNumberOfCells() << " trajectories clustered " <<std::endl;

		  if (RefinedTrajectories->GetNumberOfCells()<1)
		  {
			  std::cerr<< "The current setting of data/parameters have resulted in zero clustered trajectories"<<std::endl;
			  return EXIT_FAILURE;
		  }

		  Posterior = ComputePosterior(Likelihood,Prior);
		  UpdateModelParameters(DissimilarityMatrix,Posterior,alpha,beta,Prior,havePrior);

		  //EM Block Ends

		  if (debug)
		  {
			  std::cout << DissimilarityMatrix.max_value() << " " << DissimilarityMatrix.min_value() << std::endl;
			  std::string tempFilename;
			  tempFilename = OutputDirectory + "/trajectories_iter" + currentIteration.str() +".vtp";
			  copyField.FA = 0; copyField.Tensor = 0; copyField.Correspondences =1; copyField.CaseName=0;
			  WriteVTKfile(RefinedTrajectories, tempFilename, copyField);
			  //std::cout<< Posterior << std::endl;
			  std::cout<< "alpha = " << alpha << std::endl;
			  std::cout<< "beta = " << beta << std::endl;

			  tempFilename = OutputDirectory + "/posterior_iteration" + currentIteration.str() +".csv";
			  WriteCSVfile(tempFilename, Posterior);
		  }

		  //Update centers:
		  CenterType NewCenters = UpdateCenters(RefinedTrajectories, Centers, Posterior, MinPost);

		  spaceResolution = std::max(tractographyStepSize,(VariableType) 1.0);  //1mm space resolution >= tractographyStepSize
		  std::cout << "Setting the space resolution to " << spaceResolution << std::endl;
		  samplesDistance = std::max((VariableType) 2.5* spaceResolution,(VariableType) 2.5); //mm

		  //Smooth centers:
		  CenterType SmoothedCenters = SmoothMeshes(NewCenters, samplesDistance, 0);

		  if(debug)
		  {
			  for (unsigned int cc=0; cc<Centers.size(); cc++)
			  {
				  copyField.FA = 0; copyField.Tensor = 0; copyField.Correspondences =0; copyField.CaseName=0;
				  std::string fn = OutputDirectory + "/temp_center.vtp";
				  WriteVTKfile(SmoothedCenters.at(cc), fn, copyField);
			  }
		  }

		  Centers = SmoothedCenters;
		  Trajectories = RefinedTrajectories;
		}

	  AssignClusterLabels(Trajectories,Posterior);
	  //Done with clustering.
	  if (!analysis)
	  {
	    copyField.FA = 0; copyField.Tensor = 1; copyField.Correspondences =1; copyField.CaseName=0; copyField.ClusterLabel = 1;
	    WriteVTKfile(Trajectories, outputClustersFilename.c_str(),copyField);
	  }
	  else
	  {
	   //Start Quantitative Analysis:
	   //Compute and add diffusion scalar measures to each point on the mesh:
		  ComputeScalarMeasures(Trajectories);

	      AddOrientation(Centers);
	      cosAngleThreshold = 0.9;

		  DissimilarityMatrix = ComputeDissimilarity(Trajectories, Centers, subSpace, MaxDist, cosAngleThreshold, 1);

		  //Write the clustering output with scalars
		  copyField.FA = 1; copyField.Tensor = 0; copyField.Correspondences =1; copyField.CaseName=0; copyField.ClusterLabel = 1;
		  WriteVTKfile(Trajectories, outputClustersFilename.c_str(),copyField);


		  //Generate separate mesh for each cluster -> cluster + center:
		  std::vector <MeshType::Pointer> cluster;
		  std::vector<unsigned long int> cellId;
		  std::vector <Array3DType> clusterFeatures;

		  std::vector<std::string> clusternames = generateFilenames(centersFilename,ncells);

		  Array3DType posts;

		  for (unsigned int k=0; k<Centers.size(); ++k)
		  {
			  std::stringstream currentClusterId;
			  currentClusterId<<k+1;

			  //Separate cluster k'th
			  cluster.push_back(getCluster(Trajectories,k));

			  //Separate posterior probabilities
			  posts.push_back(getClusterPosterior(Posterior,Trajectories,k));

			  //Compute the feature matrices
			  clusterFeatures.push_back(BuildFeatureMatrix(cluster[k],Centers.at(k),k));
			  if (clusterFeatures[k].at(0).rows()>0) //not empty
			  {
				  //Now compute the mean FA and assign it to the pointvalue.FA of each point on the center
				  ArrayType meanFA; //, stdFA;
				  std::vector<unsigned int> nids;
				  meanFA = meanMat(clusterFeatures[k].at(0),-1);                          //TODo: compute the weighted average
				  //Add point data to the cell with the cell ID of cellId in the oldCenters mesh to get visualized with the
				  // mean FA after loading in Slicer:
				  AddPointScalarToACell(Centers.at(k),0,meanFA );//Centers gets updated.
			  }
			  if (!outputCentersFilename.empty())
			  {
				  copyField.FA = 1; copyField.Tensor = 0; copyField.Correspondences = 0; copyField.CaseName=0;
				  WriteVTKfile(Centers, outputCentersFilename ,copyField);
			  }

			  //write out the posterior for each cluster
			  std::string filename;
			  if (posts[k].rows()>0)
			  {
				  filename = OutputDirectory+"/cluster"+ currentClusterId.str() + "_posterior.csv";
				  WriteCSVfile(filename, posts[k]);
				  //now write out the feature matrices to files:
				  filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_FA.csv";
				  WriteCSVfile(filename,clusterFeatures[k].at(0));

				  filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_MD.csv";
				  WriteCSVfile(filename,clusterFeatures[k].at(1));

				  filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_ParDiff.csv";
				  WriteCSVfile(filename,clusterFeatures[k].at(2));

				  filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_PerDiff.csv";
				  WriteCSVfile(filename,clusterFeatures[k].at(3));

				  // Write individual files for each cluster and its center.
				  filename = OutputDirectory+ "/"+ clusternames.at(k) + "_center" + currentClusterId.str()+".vtp";
				  copyField.FA = 1; copyField.Tensor = 0; copyField.Correspondences = 0; copyField.CaseName=0;
				  WriteVTKfile(Centers.at(k),filename,copyField);

				  filename = OutputDirectory+"/"+ clusternames.at(k) + "_cluster" + currentClusterId.str()+".vtp";
				  copyField.FA = 0; copyField.Tensor = 1; copyField.Correspondences = 1; copyField.CaseName=0;
				  WriteVTKfile(cluster[k], filename,copyField);
			  }

			  std::vector<std::string> subjectNames;
			  std::vector<MeshType::Pointer> subClusters;

			  if (population)

			  {
				  std::vector<std::string> subjectNamesInCluster =  getClusterSubjectNames(cluster[k]);
				  std::string myfilename = OutputDirectory + "/cluster" + currentClusterId.str()+ "_subjectNames.txt";
				  ofstream myfile;
				  std::cout<< "Writing out "<< myfilename.c_str() << "..." << std::endl;
				  myfile.open (myfilename.c_str());
				  for (unsigned int m =0; m<subjectNamesInCluster.size(); m++)
				  {
					  myfile << subjectNamesInCluster.at(m);
					  myfile << std::endl;
				  }
				  myfile.close();
			      std::vector<Array3DType> clusterSubjectFeatures;
			      MeshType::Pointer subjectCenter = Centers.at(k);

				  for (unsigned int sn=0; sn<allfilenames.size(); sn++)
				  {
					  std::string filename = itksys::SystemTools::GetFilenameName(allfilenames.at(sn));
					  std::string subjectName = filename.substr(0,13);
					  std::string outputfilename = OutputDirectory + "/cluster" + currentClusterId.str() + "_" + subjectName + ".vtp";

					  //Generate separate meshes for each subject in the population
					  subClusters.push_back(getCluster(cluster[k], allfilenames.at(sn)));
					  copyField.FA = 0; copyField.Tensor = 1; copyField.Correspondences = 1; copyField.CaseName=0;
					  WriteVTKfile(subClusters[sn], outputfilename.c_str() ,copyField);
	       			  clusterSubjectFeatures.push_back(BuildFeatureMatrix(subClusters[sn],Centers.at(k),k));

	       			  Array2DType clusterSubjectFeatureMeans;
	       			  clusterSubjectFeatureMeans.set_size(4,subjectCenter->GetNumberOfPoints());
	       			  ArrayType meanFA = meanMat(clusterSubjectFeatures[sn].at(0),-1);
	       			  clusterSubjectFeatureMeans.set_row(0,meanMat(clusterSubjectFeatures[sn].at(0),-1));
	       			  clusterSubjectFeatureMeans.set_row(1,meanMat(clusterSubjectFeatures[sn].at(1),-1));
	       			  clusterSubjectFeatureMeans.set_row(2,meanMat(clusterSubjectFeatures[sn].at(2),-1));
	       			  clusterSubjectFeatureMeans.set_row(3,meanMat(clusterSubjectFeatures[sn].at(3),-1));

	    	  		  filename = OutputDirectory+"/cluster" + currentClusterId.str()+ "_" + subjectName +".csv";
	    	  		  WriteCSVfile(filename,clusterSubjectFeatureMeans);

	    	  		  AddPointScalarToACell(subjectCenter,0, meanFA);//Centers gets updated.
	    	  		  filename = OutputDirectory+"/center" + currentClusterId.str()+ "_" + subjectName +".vtp";
	    	  	      copyField.FA = 1; copyField.Tensor = 0; copyField.Correspondences = 0; copyField.CaseName=0;
	    	  		  WriteVTKfile(subjectCenter,filename,copyField);

				  }
			  }//done if population
			  }//done for each cluster
		  }//done if analysis
	  std::cout << "Done." << std::endl;
	  return EXIT_SUCCESS;
	}
	else //zero number of trajectories or centers
	{
		std::cout << "Completed with Error." << std::endl;
		return EXIT_FAILURE;
	}
}
