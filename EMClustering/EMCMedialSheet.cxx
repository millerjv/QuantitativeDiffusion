#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
//// SLICER VERSION
#include "Common.h"
#include "EMClusteringIO.h"
#include "Registration.h"
#include "myMaths.h"
#include "EMCore.h"
#include "MeshOperations.h"
#include "Quantification.h"

#include "EMClusteringCLP.h"

int main(int argc, char* argv[])
{
	PARSE_ARGS;

	MeshType::Pointer    Trajectories, Centers, atlasCenters;
	MeshType::Pointer    oldCenters = MeshType::New();
	ImageType::Pointer subSpace;

	//population=0;
	//use_atlas = 0;
	//analysis = 1;
	//
	bool debug = 0;
	CopyFieldType copyField = {0,0,1,1,0,1};

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

	// Get the intialcenter(s)
	if (use_atlas)
	{
		std::string path = atlas_directory;
		if (atlas_directory[0] == '.')
		{
			path = itksys::SystemTools::GetFilenamePath(argv[0]);
			path = path + "/" + atlas_directory;
		}
		std::string atlasFilename;
		atlasFilename = path + "/atlasCenters.vtp";
		atlasCenters = ReadVTKfile(atlasFilename);

		//Read FA map of the atlas
		atlasFilename = path +"/atlasFAMap.nhdr";
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		std::cout << "Reading " << atlasFilename << " ..." <<std::endl;
		reader->SetFileName(atlasFilename);
		reader->Update();
		ImageType::Pointer atlasFAVolume = reader->GetOutput();
		//Read FA map of the case
		ReaderType::Pointer reader2 = ReaderType::New();
		reader2->SetFileName(subjectFAFilename.c_str());
		std::cout << "Reading " << subjectFAFilename.c_str() << " ..." <<std::endl;
		reader2->Update();
		ImageType::Pointer caseFAVolume = reader2->GetOutput();
		//Read which bundles are selected to be clustered:           /// needs to be standardized
		std::vector<unsigned long int> atlasCellIDs;
		if (genu)
		{atlasCellIDs.push_back(0);}//{atlasCellIDs.push_back(2);}
		if (splenium)
		{atlasCellIDs.push_back(1);}//{atlasCellIDs.push_back(0);}
		if (cingulumR)
		{atlasCellIDs.push_back(3);}
		if (cingulumL)
		{atlasCellIDs.push_back(4);}

		if (atlasCellIDs.size()>0)
		{
			//Do affine registration
			TransformType::Pointer transform;
			std::cout << "Registering the atlas' FA volume to subject's ..." << std::endl;
			transform = doAffineRegistration(caseFAVolume, atlasFAVolume, OutputDirectory);
			//Select and transfer the centers
			Centers = applyTransform(atlasCenters, transform, atlasCellIDs );
			CopyFieldType copyField = {0,0,0,0};
			WriteVTKfile(Centers,transformedCentersFilename.c_str(),copyField);
		}
		else
		{
			std::cerr << "User should select at least one of the bundles in the atlas' list to be clustered! " << std::endl;
			return -1;
		}
	}
	else if (!centersFilename.empty())
	{
		Centers = ReadVTKfile(centersFilename.c_str());
	}
	else
	{
		std::cerr << "No initial center was given" << std::endl;
		return -1;
	}

	//////////////////////////////////////////////////////
	// Continue if you have valid trajectories and centers:
	//////////////////////////////////////////////////////
	if (Centers->GetNumberOfCells()>0 && Trajectories->GetNumberOfCells()>0)
	{
		VariableType tractographyStepSize  = getSampleSpacing(Trajectories);
		std::cout << "Calculated tractography step size is " << tractographyStepSize << std::endl;

		// set the space to the limits of input trajectories
	    VariableType spaceResolution = 4; //mm  initial value for the first iteration -- space resolution >= tractographyStepSize
	    std::cout << "Setting the initial space resolution to " << spaceResolution << std::endl;
	    samplesDistance = 5; //mm  initial value for the first iteration
	    // Resample initial centers:
	    bool justResample =1;
	    Centers = SmoothMesh(Centers, samplesDistance, justResample);


		VariableType MinPost = (VariableType) 1/(Centers->GetNumberOfCells());
		// If the number of clusters is 1:
		if (MinPost>0.5) { MinPost = 0.5;}

		VariableType MinLike = 0.1*MinLikelihoodThr;

		//EM Initialization

		ArrayType alpha, beta, MyMinLikelihoodThr,dd;
		Array2DType DissimilarityMatrix, Likelihood, Prior, Posterior; //NxK

		alpha.SetSize(Centers->GetNumberOfCells()); alpha.fill(1);
		beta.SetSize(Centers->GetNumberOfCells()); beta.fill(10);
		Prior.SetSize(Trajectories->GetNumberOfCells(),Centers->GetNumberOfCells());
		bool havePrior = 0;
		if (havePrior)
		{
			setPriorInfo(Prior, Trajectories);
		}
		else //in absence of an atlas
		{
			VariableType initp = 1.0 / (Centers->GetNumberOfCells());
			Prior.fill(initp);
		}

		for (int i=0; i<maxNumberOfIterations; ++i)
		{
		  std::cout<< "Iteration  " << i+1 << std::endl;
		  std::stringstream currentIteration;
		  currentIteration << i+1;

		  subSpace = getSubSpace(Trajectories, spaceResolution);

		  DissimilarityMatrix = ComputeDissimilarity(Trajectories, Centers, subSpace, samplesDistance, tractographyStepSize, MaxDist);

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
		  MeshType::Pointer NewCenters = UpdateCenters(RefinedTrajectories, Centers, Posterior, MinPost);

		  //Smooth centers:
		  MeshType::Pointer SmoothedCenters = SmoothMesh(NewCenters, samplesDistance, 0);

		  if (debug)
		  {
			  std::string tempFilename;
			  tempFilename = OutputDirectory + "/centers_iteration" + currentIteration.str() +".vtp";
			  copyField.FA = 0; copyField.Tensor = 0; copyField.Correspondences =0; copyField.CaseName=0;
			  WriteVTKfile(NewCenters,tempFilename,copyField);
			  tempFilename = OutputDirectory + "/smoothed_centers_iteration" + currentIteration.str() +".vtp";
			  copyField.FA = 0; copyField.Tensor = 0; copyField.Correspondences =0; copyField.CaseName=0;
			  WriteVTKfile(SmoothedCenters,tempFilename,copyField);
		  }

		  oldCenters = Centers;
		  Centers = SmoothedCenters;
		  Trajectories = RefinedTrajectories;
		  if (i>1)
		  {
			  dd = diffMeshes(oldCenters, Centers);
			  //std::cout<< "Difference between new centers and old ones is "<< dd.max_value() <<std::endl;
			  if (dd.max_value()<5) break;
		  }
		  spaceResolution = std::max(tractographyStepSize,(VariableType) 1.0);  //1mm space resolution >= tractographyStepSize
		  std::cout << "Setting the space resolution to " << spaceResolution << std::endl;
		  samplesDistance = std::max((VariableType) 2.5* spaceResolution,(VariableType) 2.5); //mm
		}

	  AssignClusterLabels(Trajectories,Posterior);
	  copyField.FA = 0; copyField.Tensor = 1; copyField.Correspondences =1; copyField.CaseName=0; copyField.ClusterLabel = 1;
	  WriteVTKfile(Trajectories, outputClustersFilename.c_str(),copyField);

	  //Done with clustering.
	  //////////////////////////////////////////////////////////////////////
	  //Start Quantitative Analysis:
	  //////////////////////////////////////////////////////////////////////
	  if (analysis)
	  {
		  //Compute and add diffusion scalar measures to each point on the mesh:
		  ComputeScalarMeasures(Trajectories);

		  //Generate separate mesh for each cluster -> cluster + center:
		  std::vector <MeshType::Pointer> cluster, center, centerWithData;
		  std::vector<unsigned long int> cellId;
		  std::vector <Array3DType> clusterFeatures;

		  Array3DType posts;

		  for (unsigned int k=0; k<Centers->GetNumberOfCells(); ++k)
		  {
			  std::stringstream currentClusterId;
			  currentClusterId<<k+1;

			  //Separate cluster k'th
			  cluster.push_back(getCluster(Trajectories,k));

			  //Separate center  k'th
			  cellId.clear(); cellId.push_back(k);
			  center.push_back(getTrajectories(oldCenters,cellId));

			  //Separate posterior probabilities
			  posts.push_back(getClusterPosterior(Posterior,Trajectories,k));

			  //Compute the feature matrices
			  clusterFeatures.push_back(BuildFeatureMatrix(cluster[k],center[k],k));
			  if (clusterFeatures[k].at(0).rows()>0) //not empty
			  {
				  //Now compute the mean FA and assign it to the pointvalue.FA of each point on the center
				  ArrayType meanFA; //, stdFA;
				  meanFA = meanMat(clusterFeatures[k].at(0),-1);                          //TO Do: compute the weighted average
				  //Add point data to the cell with the cell ID of cellId in the oldCenters mesh to get visualized with the
				  // mean FA after loading in Slicer:
				  AddPointScalarToACell(oldCenters,k, meanFA );//oldCenters gets updated.
			  }

			  centerWithData.push_back(getTrajectories(oldCenters,cellId));

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

				  filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_PerDiff.csv";
				  WriteCSVfile(filename,clusterFeatures[k].at(2));

				  filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_ParDiff.csv";
				  WriteCSVfile(filename,clusterFeatures[k].at(3));

				  // Write individual files for each cluster and its center.
				  filename = OutputDirectory+"/center" + currentClusterId.str()+".vtp";
				  copyField.FA = 1; copyField.Tensor = 0; copyField.Correspondences = 0; copyField.CaseName=0;
				  WriteVTKfile(centerWithData[k],filename,copyField);

				  filename = OutputDirectory+"/cluster" + currentClusterId.str()+".vtp";
				  copyField.FA = 1; copyField.Tensor = 1; copyField.Correspondences = 1; copyField.CaseName=0;
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

				  for (unsigned int sn=0; sn<allfilenames.size(); sn++)
				  {
					  std::string filename = itksys::SystemTools::GetFilenameName(allfilenames.at(sn));
					  std::string subjectName = filename.substr(0,13);
					  std::string outputfilename = OutputDirectory + "/cluster" + currentClusterId.str() + "_" + subjectName + ".vtp";

					  //Generate separate meshes for each subject in the population
					  subClusters.push_back(getCluster(cluster[k], allfilenames.at(sn)));
					  copyField.FA = 1; copyField.Tensor = 1; copyField.Correspondences = 1; copyField.CaseName=0;
					  WriteVTKfile(subClusters[sn], outputfilename.c_str() ,copyField);
				  }
			  }
		  }
		  copyField.FA = 1;
	  }
	  else
	  {
	  copyField.FA = 0;
	  }
	  //write centers with new point data from the quantitative analysis part.
	  copyField.Tensor = 0; copyField.Correspondences = 0;copyField.CaseName=0;
	  WriteVTKfile(oldCenters, outputCentersFilename.c_str(),copyField);

	  std::cout << "Done." << std::endl;
	  return EXIT_SUCCESS;
	}
	else
	{
		return EXIT_FAILURE;
	}
}
