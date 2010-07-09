#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
//// SLICER VERSION
#include "EMClusteringIO.h"
#include "Registration.h"
#include "myMaths.h"

#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkBoundingBox.h>

#include "EMClusteringCLP.h"


Array2DType ComputeDissimilarity(MeshType* mesh, MeshType* mesh_centers, ImageType* space, VariableType deltaS, VariableType deltaT, VariableType maxDistance)
{
	unsigned long int NumberOfTrajectories=mesh->GetNumberOfCells();
	unsigned int NumberOfClusters=mesh_centers->GetNumberOfCells();
	Array2DType DissimilarityMeasure;
	DissimilarityMeasure.SetSize(NumberOfTrajectories,NumberOfClusters);
	VariableType LargeDist = itk::NumericTraits<VariableType>::max();
	DissimilarityMeasure.fill(LargeDist);
	VariableType maxDist = maxDistance/space->GetSpacing()[0];  //Conversion from mm to space's units
	bool considerOrientation = 0;

	for (unsigned int ClusterIdx = 0; ClusterIdx<NumberOfClusters; ++ClusterIdx)
	{
		// Build a distance map for each cluster center
		space->FillBuffer(0);
		int myCurrentLabel =0;

		// Mark up the voxels of the centerlines on the image
		CellAutoPointer Centerline;
		mesh_centers->GetCell(ClusterIdx, Centerline);
		PolylineType::PointIdIterator pit = Centerline->PointIdsBegin();
		ImageType::IndexType ind;
		MeshType::PointType point, nextPoint, prevPoint;
		std::vector<long int> MyLabels;
		std::vector <VectorType> orientationOnCenter;
		MyLabels.clear();
		bool outOfSpace = 1;
		for (unsigned int j=0; j < Centerline->GetNumberOfPoints(); ++j)
		{
			mesh_centers->GetPoint(*pit, &point);
			if (space->TransformPhysicalPointToIndex(point, ind))
			{
				outOfSpace = 0;
				myCurrentLabel++;
				space->SetPixel(ind,myCurrentLabel);

				//populate the labels in MyLabels along each center
				MyLabels.push_back(myCurrentLabel);

				if (considerOrientation)
				{
					if (j==0)
					{
						mesh_centers->GetPoint(*(pit+1), &nextPoint);
						prevPoint = point;
					}
					else if (j==Centerline->GetNumberOfPoints()-1)
					{
						mesh_centers->GetPoint(*(pit-1), &prevPoint);
						nextPoint = point;
					}
					else
					{
						mesh_centers->GetPoint(*(pit-1), &prevPoint);
						mesh_centers->GetPoint(*(pit+1), &nextPoint);
					}
					VectorType orientationPoint = (nextPoint - prevPoint);
					orientationPoint.Normalize();
					orientationOnCenter.push_back(orientationPoint);
				}
			}
			else
			{
				std::cout<<"Point "<< point<<" on the center "<< ClusterIdx+1  <<" is out of the space"<<std::endl;
			}
			++pit;
		}

		if (!outOfSpace)
		{
			std::cout<< " Generating the Distance Map for Cluster "<< ClusterIdx+1 <<" ..."<< std::endl;

			// Apply the daneilsson filter
			typedef itk::DanielssonDistanceMapImageFilter< ImageType, ImageType > DistanceMapFilterType;
			DistanceMapFilterType::Pointer DMFilter = DistanceMapFilterType::New();
			DMFilter->SetInput(space);
			DMFilter->InputIsBinaryOff();
			DMFilter->Update();

			ImageType::Pointer DistanceMap=DMFilter->GetOutput();
			ImageType::Pointer VoronoiMap=DMFilter->GetVoronoiMap();

			//write out the images:
			/*typedef itk::ImageFileWriter< ImageType > WriterType;
      	  	  WriterType::Pointer writer = WriterType::New();
      	  	  writer->SetFileName("myspace.nhdr");
      	  	  writer->SetInput(VoronoiMap);
      	  	  std::cout<< "Writing the voronoi map" << std::endl;
      	  	  writer->Update(); */

			//Create the interpolator for the Distance Map
			itk::LinearInterpolateImageFunction<ImageType, CoordinateType>::Pointer DMInterpolator =
					itk::LinearInterpolateImageFunction<ImageType, CoordinateType>::New();
			DMInterpolator->SetInputImage(DistanceMap);

			//Find the dissimilarity measure for each trajectory
			for (unsigned int t=0; t<NumberOfTrajectories; ++t)
			{
				CellAutoPointer atrajectory;
				mesh->GetCell(t, atrajectory);
				PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
				ImageType::IndexType ind;
				MeshType::PointType point, nextPoint, prevPoint;
				MeshType::PixelType pointvalue;
				VariableType sumdist = 0, sumSinAngle=0;
				VariableType currentPointDist;
				std::vector<long int> MyLabelsOnTrajectory;
				std::vector<VectorType> orientationOnTrajectory;
				orientationOnTrajectory.clear();
				MyLabelsOnTrajectory.clear();
				long int lastLabel = 0;
				unsigned int NumberOfPointsOnTrajectory = atrajectory->GetNumberOfPoints();

				for (unsigned int j=0; j < NumberOfPointsOnTrajectory; ++j)
				{
					mesh->GetPoint(*pit, &point);
					if (space->TransformPhysicalPointToIndex(point, ind))
					{
						itk::ContinuousIndex<double, 3> cidx;
						if (space->TransformPhysicalPointToContinuousIndex( point, cidx ))
						{
							currentPointDist = DMInterpolator->Evaluate(point); //Note : not in mm
							sumdist+= currentPointDist;
						}
						mesh->GetPointData( *pit, &pointvalue );
						if (considerOrientation)
						{
							if (j==0)
							{
								mesh_centers->GetPoint(*(pit+1), &nextPoint);
								prevPoint = point;
							}
							else if (j==NumberOfPointsOnTrajectory-1)
							{
								mesh_centers->GetPoint(*(pit-1), &prevPoint);
								nextPoint = point;
							}
							else
							{
								mesh_centers->GetPoint(*(pit-1), &prevPoint);
								mesh_centers->GetPoint(*(pit+1), &nextPoint);
							}
							VectorType orientationPoint = (nextPoint - prevPoint);
							orientationPoint.Normalize();
							orientationOnTrajectory.push_back (orientationPoint);
						}
						pointvalue.Correspondence.set_size(NumberOfClusters);
						//1st output -- the correspondence info is going to be used in updating
						//the centers and further quantitative analysis.

						int  tempLabel = VoronoiMap->GetPixel(ind);
						if (considerOrientation)
						{
							VariableType cosAngle = orientationOnTrajectory.at(j)*orientationOnCenter.at(tempLabel);
							std::cout <<cosAngle << " " << std::endl;
							sumSinAngle += sqrt(1-cosAngle*cosAngle);
						}

						if (currentPointDist >maxDist)//(tempLabel< lastLabel || currentPointDist >maxDist)
						{
							pointvalue.Correspondence[ClusterIdx] = -100; //NaN
						}
						else
						{
							if (tempLabel != lastLabel)
							{
								MyLabelsOnTrajectory.push_back(tempLabel);
							}
							pointvalue.Correspondence[ClusterIdx] = tempLabel;
							lastLabel = tempLabel;
						}

						mesh->SetPointData( *pit, pointvalue );

					}
					else
					{
						std::cout << "Point " << point <<  " is outside the space" <<std::endl;
					}

					++pit;
				}

				VariableType AveDist = sumdist*space->GetSpacing()[0]/NumberOfPointsOnTrajectory;
				VariableType AveSinAngle = sumSinAngle/NumberOfPointsOnTrajectory;

				unsigned int NumberOfMissingPoints = MyLabels.size() - MyLabelsOnTrajectory.size();
				//2nd output
				DissimilarityMeasure[t][ClusterIdx] = AveSinAngle + AveDist*(1+ (NumberOfMissingPoints*deltaS)/(NumberOfPointsOnTrajectory*deltaT));
				//std::cout << " " << sumdist << " " << atrajectory->GetNumberOfPoints() << " "<< AveDist<<" "<< MyLabels.size()<< " " <<(NumberOfMissingPoints*deltaS)/(NumberOfPointsOnTrajectory*deltaT)<<" " << DissimilarityMeasure[t][ClusterIdx]<<std::endl;
			}
		}
		else
		{
			std::cout<< "Center "<< ClusterIdx+1 <<" is out of space."<< std::endl;
		}

	}

	return DissimilarityMeasure; //NxK
}

MeshType::Pointer RefineData(const MeshType* mesh, Array2DType &DissimilarityMatrix, Array2DType &Likelihood, Array2DType &Prior, ArrayType MyMinLikelihoodThr, bool havePrior)
{

  MeshType::Pointer RefinedMesh = MeshType::New();
  unsigned long pitNew = 0;
  unsigned long citNew = 0;
  ArrayType arow, peaks;
  Array2DType newDissimilarityMatrix, newLikelihood, newPrior;

  newDissimilarityMatrix.SetSize(DissimilarityMatrix.rows(),DissimilarityMatrix.cols());
  newLikelihood.SetSize(Likelihood.rows(),Likelihood.cols());
  newPrior.SetSize(Prior.rows(),Prior.cols());
  //Store the peak value of likelihood function for each cluster
  peaks.SetSize(Likelihood.cols());
  for (unsigned long int cc=0; cc<Likelihood.cols(); ++cc)
  {
	  peaks[cc] = Likelihood.get_column(cc).max_value();
	  //std::cout << peaks[cc] << std::endl;
  }

  for (unsigned long int n=0; n<Likelihood.rows(); ++n) //go over the trajectories
  {
    arow = Likelihood.get_row(n);
    //To Do: set the threshold differently when havePrior = 1
    bool copycell = 0;
    for (unsigned int m=0; m < arow.Size(); ++m)  //go over the clusters
    if (arow(m)> MyMinLikelihoodThr(m)||DissimilarityMatrix.get_row(n)(m)<= 5 )   //copy the trajectories with d<5mm and likelihood> minLikelihood
      { // To Do: || arow(m)<peaks(m)"
        copycell = 1;
        break;
      }

      if (copycell)
      {
        newDissimilarityMatrix.set_row(citNew, DissimilarityMatrix.get_row(n));
        newLikelihood.set_row(citNew, Likelihood.get_row(n));
        newPrior.set_row(citNew, Prior.get_row(n));


        CellAutoPointer atrajectory, copiedtrajectory;
        copiedtrajectory.TakeOwnership( new PolylineType );

        mesh->GetCell(n, atrajectory);
        PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
        MeshType::PointType point;
        MeshType::PixelType pointvalue;
        CellDataType cellvalue;

        for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
        {
          mesh->GetPoint(*pit, &point);
          mesh->GetPointData( *pit, &pointvalue );
          RefinedMesh->SetPoint(pitNew, point);
          RefinedMesh->SetPointData( pitNew, pointvalue );
          copiedtrajectory->SetPointId(j,pitNew);
          ++pit; ++pitNew;
        }

        mesh->GetCellData(n, &cellvalue);
        RefinedMesh->SetCell(citNew, copiedtrajectory);
        RefinedMesh->SetCellData(citNew, cellvalue);
        citNew++;
      }
  }

  newDissimilarityMatrix = newDissimilarityMatrix.extract(citNew, newDissimilarityMatrix.cols());
  newLikelihood = newLikelihood.extract( citNew, newLikelihood.cols());
  newPrior = newPrior.extract( citNew, newPrior.cols());

  DissimilarityMatrix.SetSize( newDissimilarityMatrix.rows(), newDissimilarityMatrix.cols());
  Likelihood.SetSize( newLikelihood.rows(), newLikelihood.cols());
  Prior.SetSize( newPrior.rows(), newPrior.cols());

  DissimilarityMatrix.copy_in(newDissimilarityMatrix.data_block());
  Likelihood.copy_in(newLikelihood.data_block());
  Prior.copy_in(newPrior.data_block());

  return RefinedMesh;
}


void UpdateModelParameters(const Array2DType &DissimilarityMatrix, const Array2DType &Posterior, ArrayType& alpha, ArrayType& beta, Array2DType& W, bool havePrior)
{

  Array2DType A,B;
  A = element_product(Posterior, DissimilarityMatrix);
  B = element_product(Posterior, DissimilarityMatrix.apply(logf));
  unsigned int NumberOfClusters = DissimilarityMatrix.cols();
  ArrayType X, N;
  X.SetSize(NumberOfClusters);
  N.SetSize(NumberOfClusters);

  for (unsigned int k=0; k<NumberOfClusters; ++k)
  {
    N[k] = Posterior.get_column(k).sum();
    if (N[k]>0)
    {
      X[k] = log(A.get_column(k).sum()/N(k))- B.get_column(k).sum()/N(k);
      alpha[k] = (3- X(k) + sqrt(pow(X(k)-3,2)+ 24*X(k)))/(12*X(k));
      beta[k] = A.get_column(k).sum()/(alpha(k)*N(k));
      if (!havePrior)
      { //update the mixing weights
        W.set_column(k, N[k]/Posterior.rows());
      }
    }
    else
    { //null cluster - no distribution
      alpha[k] = 0;
      beta[k] = 0;
      const VariableType z = 0;
      W.set_column(k,z);
    }
  }
  // std::cout << W <<std::endl;
}
void AssignClusterLabels(MeshType* mesh, const Array2DType &Posterior)
{
  if (Posterior.rows()!= mesh->GetNumberOfCells())
  {
    std::cerr<< "There is a miss-match between the number of trajectories and the membership probabilities to be assigned" <<std::endl;
  }
  CellDataType cellvalue;
  ArrayType arow;
  for (unsigned int t=0; t<Posterior.rows(); ++t)
  {
    VariableType my_max = -1;
    int my_max_idx;
    arow = Posterior.get_row(t);
    // find the maximum membership probability for t'th trajectory
    for (unsigned int m=0; m<arow.Size(); ++m)
    {
      if (arow(m)>my_max)
      {
        my_max = arow(m);
        my_max_idx = m;
      }
    }
    CellAutoPointer atrajectory;
    mesh->GetCell(t, atrajectory);
    mesh->GetCellData(t, &cellvalue );
    cellvalue.ClusterLabel = my_max_idx+1;    //start the Cluster Labels (Ids) from 1
    cellvalue.membershipProbability.SetSize(arow.Size());
    cellvalue.membershipProbability = arow;
    mesh->SetCellData(t, cellvalue );
  }
}

MeshType::Pointer UpdateCenters(MeshType* mesh, MeshType* mesh_centers, const Array2DType &Posterior, VariableType MinPosterior)
{
  MeshType::Pointer mesh_newcenters=MeshType::New();
  ArrayType post;
  unsigned int NumberOfClusters=mesh_centers->GetNumberOfCells();
  long int cit =0;

  for (unsigned int k=0; k<NumberOfClusters; ++k)
  {
    post = Posterior.get_column(k);
    std::vector<long int> sigIDs;
    for (unsigned long int p=0; p<post.Size(); ++p)
    {
      if (post(p)>MinPosterior)
      {
        sigIDs.push_back(p);
      }
    }

    if(sigIDs.size()<1)
    {
    	std::cout<< "Center " << k+1 << " is not being updated" <<std::endl;
    }

    CellAutoPointer Centerline, new_center;
    new_center.TakeOwnership( new PolylineType );
    mesh_centers->GetCell(k,Centerline);
    PolylineType::PointIdIterator mcit = Centerline->PointIdsBegin();
   MeshType::PointType point, last_mean_point;
    VariableType minDistBetweenSuccessivePointsOnCenter =0.5;    //the samples on the centers do not to be closer than 0.5mm
    VariableType distBetweenSuccessivePoints;

    int MyLabel, currentLabel=0;
    unsigned int s=0;
    for (unsigned int c=0; c<Centerline->GetNumberOfPoints(); ++c)
    {
      mesh_centers->GetPoint(*mcit, &point);
      /*if(refImage->TransformPhysicalPointToIndex(point, ind))
      {*/
      MyLabel = ++currentLabel;
      MeshType::PixelType tpointvalue;
   	  std::vector<MeshType::PointType> pntStack;
   	  std::vector<VariableType> postStack;

   	  MeshType::PointType tpoint, mean_point,sum_points, closest_point;
      VariableType dist,closest_point_post;
      pntStack.clear();
      postStack.clear();
   	  if (sigIDs.size()>0)
   	  {//update the center by taking the average of trajectories
    	for (unsigned long int t=0; t<sigIDs.size(); ++t)
    	{
    	  CellAutoPointer atrajectory;
    	  mesh->GetCell(sigIDs.at(t),atrajectory);
		  PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
		  VariableType MinDist = itk::NumericTraits<VariableType>::max();
    	  for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
    	  {
    		mesh->GetPoint(*pit, &tpoint);
    		mesh->GetPointData( *pit, &tpointvalue );
    		 //find the points on a single trajectory that corresponds to the current point on the center
    		if (tpointvalue.Correspondence(k)==MyLabel)
    		{
    		  dist = tpoint.EuclideanDistanceTo(point);
    		  if (dist<MinDist)
    		  {
    			  MinDist = dist;
    			  closest_point = tpoint;
    			  closest_point_post = post(sigIDs.at(t));
    		  }
    		}
    		pit++;
    	  }//for j
    	  pntStack.push_back(closest_point);
          postStack.push_back(closest_point_post);

    	}//for t
        sum_points.Fill(0);
    	VariableType SumPost = 0;
    	for ( unsigned int m=0; m<pntStack.size(); ++m)
    	{
    	  MeshType::PointType temp = pntStack[m];
    	  sum_points[0] += temp[0]*postStack[m];
    	  sum_points[1] += temp[1]*postStack[m];
    	  sum_points[2] += temp[2]*postStack[m];
    	  SumPost +=postStack[m];
    	}
        if (SumPost>0)
        {
        	mean_point[0] = sum_points[0]/SumPost;
    		mean_point[1] = sum_points[1]/SumPost;
    		mean_point[2] = sum_points[2]/SumPost;
			//compute the distance between the current mean point and the previous one:
    	    if (c>0 && cit>0) //not at the beginning of the centerline
    	    {
    	    	distBetweenSuccessivePoints = mean_point.EuclideanDistanceTo(last_mean_point);
    		}
    	    else
    	    {
    		    distBetweenSuccessivePoints = minDistBetweenSuccessivePointsOnCenter;
    	    }

			if (distBetweenSuccessivePoints>= minDistBetweenSuccessivePointsOnCenter)
    		{
    		    mesh_newcenters->SetPoint(cit,mean_point);
    			new_center->SetPointId(s,cit);
    			last_mean_point = mean_point;
    			++cit; ++s;
    		}
        }//if Sumpost<0
        else //NO MATCHING POINT
    	{
    	  std::cout<<"A point on the center is not being updated!" << std::endl;
    	}
    	++mcit;
   	  }
   	  else //NO TRAJECTORIES
   	  {//just copy the point from the center to newcenters

   		  mesh_newcenters->SetPoint(cit,point);
    	  new_center->SetPointId(s,cit);
    	  ++cit; ++s; ++mcit;
      }
    }//end for c -- point along the center
    mesh_newcenters->SetCell(k,new_center);
  }//end for k (cluster)
  return mesh_newcenters;
}


MeshType::Pointer SmoothMesh(MeshType* mesh, VariableType distanceBetweenSamples, bool justResample =0)
{
  MeshType::Pointer smoothedMesh = MeshType::New();
  unsigned int NumberOfCells = mesh->GetNumberOfCells();
  long unsigned int newpid = 0;
  for (unsigned int k=0; k<NumberOfCells; ++k)
  {
    CellAutoPointer aCell,newCell;
    newCell.TakeOwnership( new PolylineType );

    mesh->GetCell(k,aCell);

    MeshType::PointType point;
    CurveType MyCurve, SmoothedCurve;
    MyCurve.SetSize(aCell->GetNumberOfPoints(),3);

    PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
    for (unsigned int j=0; j < aCell->GetNumberOfPoints(); ++j)
    {
      mesh->GetPoint(*pit, &point);
      MyCurve.set_row(j, point.GetVnlVector());
      pit++;
    }
    if (justResample)
    {
    	SmoothedCurve = SmoothAndResampleCurve(MyCurve, distanceBetweenSamples);
    }
    else
    {
    	SmoothedCurve = SmoothCurve(MyCurve,distanceBetweenSamples);
    	SmoothedCurve = SmoothAndResampleCurve(SmoothedCurve, distanceBetweenSamples);
    }
    //Put the new curve in the mesh:
    for (unsigned int j=0; j < SmoothedCurve.rows(); ++j)
    {
      MeshType::PointType mpoint;
      mpoint[0] = SmoothedCurve(j,0);
      mpoint[1] = SmoothedCurve(j,1);
      mpoint[2] = SmoothedCurve(j,2);
      smoothedMesh->SetPoint(newpid, mpoint);
      newCell->SetPointId(j,newpid);
      newpid++;
    }
    smoothedMesh->SetCell(k,newCell);
  }
  return smoothedMesh;
}

ArrayType diffMeshes(const MeshType* mesh1, const MeshType* mesh2)
{
  ArrayType dist;

  unsigned int NumberOfCells1 = mesh1->GetNumberOfCells();
  unsigned int NumberOfCells2 = mesh2->GetNumberOfCells();
  dist.SetSize(NumberOfCells1); dist.fill(10000);
  if (NumberOfCells1 != NumberOfCells2)
  {
    std::cout << "Number of cells don't match between the given meshes" <<std::endl;
  }
  else
  {
    dist.fill(0);
    for (unsigned int k=0; k<NumberOfCells1; ++k)
    {
      CellAutoPointer aCell1,aCell2;
      //newCell.TakeOwnership( new PolylineType );

      mesh1->GetCell(k,aCell1);
      mesh2->GetCell(k,aCell2);

      MeshType::PointType point;
      CurveType MyCurve1, MyCurve2;
      MyCurve1.SetSize(aCell1->GetNumberOfPoints(),3);
      MyCurve2.SetSize(aCell2->GetNumberOfPoints(),3);

      PolylineType::PointIdIterator pit1 = aCell1->PointIdsBegin();
      for (unsigned int j=0; j < aCell1->GetNumberOfPoints(); ++j)
      {
        mesh1->GetPoint(*pit1, &point);
        MyCurve1.set_row(j, point.GetVnlVector());
        pit1++;
      }

      PolylineType::PointIdIterator pit2 = aCell2->PointIdsBegin();
      for (unsigned int i=0; i < aCell2->GetNumberOfPoints(); ++i)
      {
        mesh2->GetPoint(*pit2, &point);
        MyCurve2.set_row(i, point.GetVnlVector());
        pit2++;
      }


      dist[k] = diffCurve(MyCurve1, MyCurve2);

    }
  }
  return dist;

}



void SetInitialValue(const Array2DType &DissimilarityMatrix, ArrayType &beta)
{
  VariableType MaxErr = 10;
  for (unsigned int k = 0; k<DissimilarityMatrix.cols(); ++k)
  {
    VariableType sumErr = 0; long int n =0;
    for (unsigned long int t = 0; t<DissimilarityMatrix.rows(); ++t)
      if (DissimilarityMatrix(t,k)<MaxErr)
      {
        sumErr+=DissimilarityMatrix(t,k);
        n++;
      }
      beta[k] = sumErr/n;
      //std::cout<< beta[k] << std::endl;
      //beta[k] = 5;
  }
}

ArrayType AdjustThreshold(VariableType factor, ArrayType alpha, ArrayType beta)
{
  ArrayType MyMinLikelihoodThr;
  MyMinLikelihoodThr.set_size(alpha.size());
  for (unsigned int k=0; k<alpha.size(); ++k)
  {
    if (beta[k]>0) //valid distribution
    {
      VariableType GammaMode = (alpha[k]>=1)?((alpha[k]-1)*beta[k]):0.1;
      MyMinLikelihoodThr[k] = factor*Gamma(GammaMode,alpha[k],beta[k]);
    }
    else
    {
      MyMinLikelihoodThr[k] = itk::NumericTraits<VariableType>::max();
    }
  }
  return MyMinLikelihoodThr;
}

void ComputeScalarMeasures(MeshType* Trajectories)
{
  for (unsigned long int t=0; t<Trajectories->GetNumberOfPoints(); ++t)
  {
    MeshType::PixelType pointvalue;
    itk::FixedArray<double, 9 > tensor9;
    TensorPixelType tensor6;

    Trajectories->GetPointData(t, &pointvalue );
    tensor9 = pointvalue.Tensor;
    tensor6[0]=tensor9[0];
    tensor6[1]=tensor9[1];
    tensor6[2]=tensor9[2];
    tensor6[3]=tensor9[4];
    tensor6[4]=tensor9[5];
    tensor6[5]=tensor9[8];
    EigenValuesArrayType eigenvals;
    tensor6.ComputeEigenValues(eigenvals);
    pointvalue.EigenValues = eigenvals;
    pointvalue.FA = tensor6.GetFractionalAnisotropy();

    if (pointvalue.FA>=1)
    {
      pointvalue.FA = 0.5; // for now
      std::cerr << " Reached at a point with FA greater than 1" <<std::endl;
    }
    Trajectories->SetPointData(t, pointvalue );

  }
}

Array3DType BuildFeatureMatrix(const MeshType* cluster, const MeshType* center, int clusterId)
{
  Array2DType fMatrix1,fMatrix2,fMatrix3,fMatrix4;  // NxS (Number of Trajectories x Number of Samples on the Center)
  VariableType nanVal = -1;
  unsigned long int numberOfTrajectories = cluster->GetNumberOfCells();
  unsigned int numberOfSamples = center->GetNumberOfPoints();
  fMatrix1.set_size(numberOfTrajectories,numberOfSamples);
  fMatrix2.set_size(numberOfTrajectories,numberOfSamples);
  fMatrix3.set_size(numberOfTrajectories,numberOfSamples);
  fMatrix4.set_size(numberOfTrajectories,numberOfSamples);

  MeshType::PointType point;
  int MyLabel, currentLabel=0;
  for (unsigned int s=0; s<numberOfSamples; ++s)
  {
    center->GetPoint(s, &point);
    MyLabel = ++currentLabel;
    //go over trajectories
    for (unsigned long int t=0; t<numberOfTrajectories; ++t)
    {
    	CellAutoPointer atrajectory;
    	cluster->GetCell(t,atrajectory);
    	PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
    	MeshType::PointType tpoint;
    	MeshType::PixelType tpointvalue;

   		double sumFeature1 = 0;
   		double sumFeature2 = 0;
   		double sumFeature3 = 0;
   		double sumFeature4 = 0;
   		int n=0;
   		for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
   		{
   			cluster->GetPoint(*pit, &tpoint);
   			cluster->GetPointData( *pit, &tpointvalue );
   			if (tpointvalue.Correspondence(clusterId)==MyLabel)
   			{
   				sumFeature1 += tpointvalue.FA;
   				sumFeature2 += (tpointvalue.EigenValues[0]+tpointvalue.EigenValues[1]+tpointvalue.EigenValues[2])/3;
   				sumFeature3 += (tpointvalue.EigenValues[0]+tpointvalue.EigenValues[1])/2;
   				sumFeature4 += tpointvalue.EigenValues[2];
   				n++;
   			}
   			pit++;
   		}
   		if (sumFeature1>0)
   		{
   			fMatrix1[t][s] = sumFeature1/n;
   			fMatrix2[t][s] = sumFeature2/n;
   			fMatrix3[t][s] = sumFeature3/n;
   			fMatrix4[t][s] = sumFeature4/n;
   		}
   		else
   		{
   			fMatrix1[t][s] = nanVal;
   			fMatrix2[t][s] = nanVal;
   			fMatrix3[t][s] = nanVal;
   			fMatrix4[t][s] = nanVal;
   		}
    }
  } //end for
  Array3DType allFeatures;
  allFeatures.push_back(fMatrix1);
  allFeatures.push_back(fMatrix2);
  allFeatures.push_back(fMatrix3);
  allFeatures.push_back(fMatrix4);

  return allFeatures;
}

MeshType::Pointer getTrajectories(MeshType* mesh, std::vector<unsigned long int> CellIDs)
{
  MeshType::Pointer selectedCells=MeshType::New();
  CellAutoPointer aCell, MyCell;
  long int myid =0;

  for (unsigned int c=0; c<CellIDs.size(); ++c)
  {
    MyCell.TakeOwnership( new PolylineType );
    MeshType::CellIdentifier CellID = CellIDs.at(c);
    mesh->GetCell(CellID, aCell);
    CellDataType cellvalue;
    mesh->GetCellData(CellID, &cellvalue);
    PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
    MeshType::PointType point;
    MeshType::PixelType pointvalue;
    for (unsigned int j=0; j < aCell->GetNumberOfPoints(); j++)
    {
      mesh->GetPoint(*pit, &point);
      mesh->GetPointData(*pit, &pointvalue);
      selectedCells->SetPoint(myid, point );
      selectedCells->SetPointData(myid, pointvalue );
      MyCell->SetPointId(j,myid);
      pit++; myid++;
    }
    selectedCells->SetCell(c, MyCell);
    selectedCells->SetCellData(c, cellvalue);

  }

  return selectedCells;
}

MeshType::Pointer getCluster(MeshType* Trajectories,int k)
{
  MeshType::Pointer cluster;
  std::vector<unsigned long int> CellIds;
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellAutoPointer atrajectory;
    CellDataType cellvalue;
    Trajectories->GetCell(t, atrajectory);
    Trajectories->GetCellData(t, &cellvalue );
    if (cellvalue.ClusterLabel == k+1)
    {
      CellIds.push_back(t);
    }
  }
  cluster = getTrajectories(Trajectories, CellIds);
  return cluster;
}

MeshType::Pointer getCluster(MeshType* Trajectories, std::string caseName)
{
  MeshType::Pointer cluster;
  std::vector<unsigned long int> CellIds;
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellAutoPointer atrajectory;
    CellDataType cellvalue;
    Trajectories->GetCell(t, atrajectory);
    Trajectories->GetCellData(t, &cellvalue );
    if (cellvalue.CaseName.compare(caseName)==0 )
    {
      CellIds.push_back(t);
    }
  }
  cluster = getTrajectories(Trajectories, CellIds);
  return cluster;
}

Array2DType getClusterPosterior(Array2DType Posterior, MeshType* Trajectories,int k)
{
  Array2DType post;
  std::vector<long int> CellIds;
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellAutoPointer atrajectory;
    CellDataType cellvalue;
    Trajectories->GetCell(t, atrajectory);
    Trajectories->GetCellData(t, &cellvalue );
    if (cellvalue.ClusterLabel == k+1)
    {
      CellIds.push_back(t);
    }
  }
  post.set_size(CellIds.size(),1);
  for (unsigned long int c=0; c<CellIds.size(); ++c)
  {
    post.set_row(c,Posterior.get(CellIds[c],k)); //Posterior.get_row(CellIds[c])
  }
  return post;
}

std::vector<std::string> getClusterSubjectNames(MeshType* Trajectories)
{
	std::vector<std::string> subjectNames;
	for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellDataType cellvalue;
    Trajectories->GetCellData(t, &cellvalue );
    subjectNames.push_back(cellvalue.CaseName);
  }

  return subjectNames;
}

VariableType getSampleSpacing(const MeshType* Trajectories)
{
	VariableType dist;

	MeshType::PointType p1, p2;
	if (Trajectories->GetPoint(0, &p1) && Trajectories->GetPoint(1, &p2))
	{
		dist = p1.EuclideanDistanceTo(p2);
	}
	else
	{
		std::cout << "Was not able to compute the tractography step size." << std::endl;
		return 0;
	}

	return dist;
}

ImageType::Pointer getSubSpace(const MeshType* Trajectories, VariableType spacing)
{
  ImageType::Pointer subSpace = ImageType::New();
  MeshType::BoundingBoxType const *dataBoundingBox = Trajectories->GetBoundingBox();
  dataBoundingBox->ComputeBoundingBox();
  MeshType::BoundingBoxType::BoundsArrayType dataBounds = dataBoundingBox->GetBounds();

  MeshType::PointType p1,p2,p;

  p1[0] = dataBounds[0];
  p1[1] = dataBounds[2];
  p1[2] = dataBounds[4];
  p2[0] = dataBounds[1];
  p2[1] = dataBounds[3];
  p2[2] = dataBounds[5];

  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  ImageType::SizeType size;
  size[0] = ceil((p2[0]-p1[0]+5)/spacing);
  size[1] = ceil((p2[1]-p1[1]+5)/spacing);
  size[2] = ceil((p2[2]-p1[2]+5)/spacing);

  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );

  subSpace->SetRegions(desiredRegion);
  subSpace->Allocate();

  ImageType::SpacingType imSpacing;
  imSpacing[0] = spacing;
  imSpacing[1] = spacing;
  imSpacing[2] = spacing;


  subSpace->SetSpacing(imSpacing);
  subSpace->SetOrigin(p1);

  return subSpace;
}

void  setPriorInfo(Array2DType &Prior, MeshType* Trajectories)
{
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellDataType cellvalue;
    Trajectories->GetCellData(t, &cellvalue);
    Prior.set_row(t, cellvalue.atlasPriors);
  }
}

void  AddPointScalarToACell(MeshType* mesh, MeshType::CellIdentifier CellID, ArrayType f)
{
  CellAutoPointer aCell;
  mesh->GetCell(CellID, aCell);
  PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
  MeshType::PixelType pointvalue;
  for (unsigned int j=0; j < aCell->GetNumberOfPoints(); j++)
  {
    mesh->GetPointData(*pit, &pointvalue);
    pointvalue.FA = f[j];
    mesh->SetPointData(*pit, pointvalue );
    pit++;
  }

}

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
		  MeshType::Pointer SmoothedCenters = SmoothMesh(NewCenters, samplesDistance);

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
