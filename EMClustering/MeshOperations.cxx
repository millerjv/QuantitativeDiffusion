
#include "MeshOperations.h"

Array2DType ComputeDissimilarity(MeshType* mesh, CenterType mesh_centers, ImageType* space, VariableType maxDistance, bool considerOrientation)
{
	unsigned long int NumberOfTrajectories=mesh->GetNumberOfCells();
	unsigned int NumberOfClusters=mesh_centers.size();

	Array2DType DissimilarityMeasure;
	DissimilarityMeasure.SetSize(NumberOfTrajectories,NumberOfClusters);
	VariableType LargeDist = itk::NumericTraits<VariableType>::max();
	DissimilarityMeasure.fill(LargeDist);
	VariableType maxDist = maxDistance/space->GetSpacing()[0];  //Conversion from mm to space's units

	for (unsigned int ClusterIdx = 0; ClusterIdx<NumberOfClusters; ++ClusterIdx)
	{
		// Build a distance map for each cluster center
		space->FillBuffer(0);
		int myCurrentLabel =0;

		// Mark up the voxels of the centerlines on the image
		ImageType::IndexType ind;
		MeshType::PointType point, nextPoint, prevPoint;
		MeshType::PixelType cpointvalue;
		std::vector <VectorType> orientationOnCenter;

		bool outOfSpace = 1;

		for (unsigned int j=0; j < mesh_centers.at(ClusterIdx)->GetNumberOfPoints(); ++j)
		{
			mesh_centers.at(ClusterIdx)->GetPoint(j, &point);
			mesh_centers.at(ClusterIdx)->GetPointData(j, &cpointvalue);
			if (space->TransformPhysicalPointToIndex(point, ind))
			{
				outOfSpace = 0;
				myCurrentLabel++;
				space->SetPixel(ind,myCurrentLabel);

				if (considerOrientation)
				{
				VectorType vec;
				vec.SetElement(0,cpointvalue.Orientation[0]);
				vec.SetElement(1,cpointvalue.Orientation[1]);
				vec.SetElement(2,cpointvalue.Orientation[2]);
				orientationOnCenter.push_back(vec);
				}
			}
			else
			{
				std::cout<<"Point "<< point<<" on the center "<< ClusterIdx+1  <<" is out of the space"<<std::endl;
			}

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

			//Create the interpolator for the Distance Map
			typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateType> DMInterpolatorType;
			DMInterpolatorType::Pointer DMInterpolator=DMInterpolatorType::New();

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
				VariableType sumdist = 0;
				VariableType currentPointDist;
				int missOrien=0, toofar=0;

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
						pointvalue.Correspondence.set_size(NumberOfClusters);
						//1st output -- the correspondence info is going to be used in updating
						//the centers and further quantitative analysis.

						int  currentLabel = VoronoiMap->GetPixel(ind);
						if (currentPointDist>maxDist)
						{
							pointvalue.Correspondence[ClusterIdx] = -100; //NaN
							toofar++;
						}
						else
						{
							pointvalue.Correspondence[ClusterIdx] = currentLabel;

							if (considerOrientation)
							{
								VariableType cosAngle = (pointvalue.Orientation*orientationOnCenter.at(currentLabel-1));
							    if (cosAngle>-0.5 && cosAngle<0.5)
							    {
							    	pointvalue.Correspondence[ClusterIdx] = -120; //Orientation missmatch
							    	missOrien++;
							    }
							}

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
                int misses = missOrien + toofar; int valids =  NumberOfPointsOnTrajectory - misses;
				//2nd output
                if (valids>misses)
                {
                	DissimilarityMeasure[t][ClusterIdx] = AveDist*(1+ ((float) missOrien)/(NumberOfPointsOnTrajectory));
                }
                else
                {
                	DissimilarityMeasure[t][ClusterIdx] = 50;
                }
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

CenterType SmoothMeshes(CenterType meshes, VariableType distanceBetweenSamples, bool justResample =0)
{
	CenterType smoothed;
	for(unsigned int k=0; k<meshes.size(); k++)
	{
		MeshType::Pointer smoothed_mesh = SmoothMesh(meshes.at(k), distanceBetweenSamples, justResample);
	    smoothed.push_back(smoothed_mesh);
	}
	return smoothed;

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

MeshType::PointType mean_w(std::vector<MeshType::PointType> p, std::vector<VariableType> w)
{
	MeshType::PointType mp, sp;
	sp.Fill(0);
	VariableType sw=0;
	for (unsigned int n=0; n<p.size(); n++)
	{
		//std::cout << p.at(n)[0] << p.at(n)[1] << p.at(n)[2] <<std::endl;;
		sp[0] = w.at(n)*p.at(n)[0] + sp[0];
		sp[1] = w.at(n)*p.at(n)[1] + sp[1];
		sp[2] = w.at(n)*p.at(n)[2] + sp[2];
		sw+=w.at(n);
	}
	//mp = sp/p.size();
    mp[0] = sp[0]/sw;
    mp[1] = sp[1]/sw;
    mp[2] = sp[2]/sw;
    //std::cout << mp[0] << mp[1] << mp[2] <<std::endl;;
    return mp;
}

/*
void UpdateCenters(const MeshType* mesh, CenterType mesh_centers, const Array2DType &Posterior, VariableType MinPosterior)
{

	unsigned int NumberOfClusters=mesh_centers.size();
	ArrayType post;
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
		else
		{
			itk::FixedArray<std::vector<MeshType::PointType>, 1000> points;
			itk::FixedArray<std::vector<VariableType>, 1000> probs;
			MeshType::PixelType tpointvalue;
			MeshType::PointType tpoint;


			for (unsigned long int t=0; t<sigIDs.size(); ++t)
			{
				CellAutoPointer atrajectory;
				mesh->GetCell(sigIDs.at(t),atrajectory);
				PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
				for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
				{
					mesh->GetPoint(*pit, &tpoint);
					mesh->GetPointData( *pit, &tpointvalue );
					points[tpointvalue.Correspondence(k)-1].push_back(tpoint);
					probs[tpointvalue.Correspondence(k)-1].push_back(post(sigIDs.at(t)));
					pit++;
				}//for j
			}//for t

			for (unsigned int c=0; c<mesh_centers.at(k)->GetNumberOfPoints(); ++c) //for over center points
			{
				//Updating point on the surface
				MeshType::PointType point, mean_point;
				if (points[c].size()>0)
				{
					mesh_centers.at(k)->GetPoint(c, &point);
					mean_point = mean_w(points[c], probs[c]);
					mesh_centers.at(k)->SetPoint(c,mean_point);
				}
				else
				{
					std::cout << "point " << c << " not being updated" <<std::endl;
				}

			}//end for c -- point on the center
		}
	}//end for k (cluster)
}
*/



CenterType UpdateCenters(const MeshType* mesh, CenterType centers, const Array2DType &Posterior, VariableType MinPosterior)
{
	//TODO: See if there is a more efficient/faster way to do this
  CenterType newcenters;
  ArrayType post;
  unsigned int NumberOfClusters = centers.size();
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
    	newcenters.push_back(centers.at(k));
    }
    else
    {
    	MeshType::Pointer updatedcenter=MeshType::New();
    	CellAutoPointer Centerline, new_center;
    	new_center.TakeOwnership( new PolylineType );
    	centers.at(k)->GetCell(0,Centerline);  //each center mesh has only one cell.
    	PolylineType::PointIdIterator mcit = Centerline->PointIdsBegin();
    	MeshType::PointType point, last_mean_point;
    	VariableType minDistBetweenSuccessivePointsOnCenter =0.5;    //the samples on the centers do not to be closer than 0.5mm
    	VariableType distBetweenSuccessivePoints;

    	int MyLabel, currentLabel=0;
    	unsigned int s=0;
    	for (unsigned int c=0; c<Centerline->GetNumberOfPoints(); ++c)
    	{
    		centers.at(k)->GetPoint(*mcit, &point);
    		MyLabel = ++currentLabel;
            MeshType::PixelType tpointvalue;
            std::vector<MeshType::PointType> pntStack;
            std::vector<VariableType> postStack;

            MeshType::PointType tpoint, mean_point,sum_points, closest_point;
            VariableType dist,closest_point_post;
            pntStack.clear();
            postStack.clear();
            //update the center by taking the average of trajectories
            for (unsigned long int t=0; t<sigIDs.size(); ++t)
            {
            	CellAutoPointer atrajectory;
            	mesh->GetCell(sigIDs.at(t),atrajectory);
            	PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
            	VariableType MinDist = itk::NumericTraits<VariableType>::max();
            	bool foundAMatch = 0;
            	for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
            	{
            		mesh->GetPoint(*pit, &tpoint);
            		mesh->GetPointData( *pit, &tpointvalue );
            		//find the points on a single trajectory that corresponds to the current point on the center
            		if (tpointvalue.Correspondence(k)==MyLabel)
            		{
            			foundAMatch =1;
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
            	if (foundAMatch)
            	{
            		pntStack.push_back(closest_point);
            		postStack.push_back(closest_point_post);
            	}
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
            		updatedcenter->SetPoint(cit,mean_point);
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
    	}//end for c -- point along the center
    	updatedcenter->SetCell(0,new_center);
    	newcenters.push_back(updatedcenter);
    }//else
  }//end for k (cluster)
  return newcenters;
}


void AddOrientation(MeshType* mesh)
{
	for (unsigned int t=0; t<mesh->GetNumberOfCells(); ++t)
	{
		CellAutoPointer atrajectory;
		mesh->GetCell(t, atrajectory);
		PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
	    MeshType::PointType point, nextPoint, prevPoint;
		MeshType::PixelType pointvalue;
		unsigned int NumberOfPointsOnTrajectory = atrajectory->GetNumberOfPoints();

		for (unsigned int j=0; j < NumberOfPointsOnTrajectory; ++j)
		{
			mesh->GetPointData( *pit, &pointvalue );
			VectorType orientationAtPoint;
		    if (j==0)
		    {
		    	mesh->GetPoint(*(pit+1), &nextPoint);
		    	prevPoint = point;
		    }
			else if (j==NumberOfPointsOnTrajectory-1)
		    {
				mesh->GetPoint(*(pit-1), &prevPoint);
				nextPoint = point;
			}
			else
			{
				mesh->GetPoint(*(pit-1), &prevPoint);
				mesh->GetPoint(*(pit+1), &nextPoint);
			}
			orientationAtPoint = (nextPoint - prevPoint);
			orientationAtPoint.Normalize();
			pointvalue.Orientation = orientationAtPoint;
			mesh->SetPointData( *pit, pointvalue );
			++pit;
		}
	}
}

void AddOrientation(CenterType centers)
{
  for (unsigned int k=0; k<centers.size(); k++)
	  AddOrientation(centers.at(k));
}
