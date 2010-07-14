
#include "EMCore.h"
#include "myMaths.h"

void  setPriorInfo(Array2DType &Prior, MeshType* Trajectories)
{
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellDataType cellvalue;
    Trajectories->GetCellData(t, &cellvalue);
    Prior.set_row(t, cellvalue.atlasPriors);
  }
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

Array2DType ComputePosterior(const Array2DType &Likelihood, const Array2DType &Prior)
{
  Array2DType Posterior;
  Posterior = element_product(Likelihood, Prior);
  //now normalize it:
  if (Likelihood.cols() == 1)
  {
    Posterior /= Posterior.max_value();
  }

  else if (Likelihood.cols() > 1)
  {
    for (unsigned int r = 0; r<Posterior.rows(); ++r)
    {
      VariableType n = Posterior.get_row(r).sum();
      if (n>0)
      {
        Posterior.scale_row(r, 1/n);
      }
    }
  }

  //std::cout<< Posterior << std::endl;
  return Posterior; //NxK
}

Array2DType ComputeLikelihood(const Array2DType &DissimilarityMatrix, ArrayType alpha, ArrayType beta)
{
  unsigned int NumberOfClusters = DissimilarityMatrix.cols();
  unsigned int NumberOfTrajectories = DissimilarityMatrix.rows();

  Array2DType Likelihood;
  Likelihood.SetSize(NumberOfTrajectories,NumberOfClusters);

  for (unsigned int k=0; k<NumberOfClusters; ++k)
  {
    if (beta[k]>0) //valid distribution
    {
      for (unsigned long int n=0; n<NumberOfTrajectories; ++n)
      {
        VariableType x = DissimilarityMatrix(n,k);
        Likelihood[n][k] = Gamma(x,alpha[k],beta[k]);
      }
    }
    else
    {
      const VariableType z = 0;
      Likelihood.set_column(k,z);
    }
  }
  return Likelihood; //NxK
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
