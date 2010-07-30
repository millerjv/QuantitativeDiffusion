
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

