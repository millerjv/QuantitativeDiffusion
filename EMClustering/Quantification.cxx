
#include "Quantification.h"

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
