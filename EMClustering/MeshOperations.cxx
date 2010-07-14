
#include "MeshOperations.h"

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
