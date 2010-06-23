#include "EMClusteringIO.h"
#include <myMaths.h>
#include <itkThinPlateSplineKernelTransform.h>

ArrayType meanMat(const Array2DType &X, int nanVal=0)
//take the column-wise mean of the matrix X, ignoring the zero elements.
{
  ArrayType mX;
  mX.SetSize(X.cols());
  ArrayType aCol;
  for (unsigned int c = 0; c<X.cols(); c++)
  {

    aCol = X.get_column(c);

    // if there is no NAN in the column:   (here we assume that the nanVal is 0 or a negative number.
    if (aCol.min_value() != nanVal)
    {
      mX(c) = aCol.mean();
    }
    else
    {
      double s = 0;
      unsigned int n = 0;
      for (unsigned int i =0; i<aCol.Size(); i++)
      {
        if (aCol(i)!= nanVal)
        {
          s+=aCol(i);
          n++;
        }
      }
      if (n!=0)
      {
        mX(c) = s/n;
      }
      else
      {
        //
        mX(c)=0;
        std::cout << "NaN column at " << c << "!" <<std::endl;
      }
    }

  }

  return mX;
}

ArrayType stdMat(const Array2DType &X, int nanVal=0)
//take the column-wise std of the matrix X, ignoring the nonVal elements.
{
  ArrayType mX;
  mX.SetSize(X.cols());
  ArrayType aCol;
  for (unsigned int c = 0; c<X.cols(); c++)
  {

    aCol = X.get_column(c);

    // if there is no NAN in the column:   (here we assume that the nanVal is 0 or a negative number.
    if (aCol.min_value() != nanVal)
    {
		mX(c) = (aCol - aCol.mean()).rms();
    }
    else
    {
      double s = 0;
      unsigned int n = 0;
      for (unsigned int i =0; i<aCol.Size(); i++)
      {
        if (aCol(i)!= nanVal)
        {
          s+=aCol(i);
          n++;
        }
      }
      if (n!=0)
      {
         double meanVal =  s/n;
		 double ms = 0;
		 for (unsigned int i =0; i<aCol.Size(); i++)
			{
			if (aCol(i)!= nanVal)
			{
				ms+=(aCol(i)-meanVal)*(aCol(i)-meanVal);

			}
			}


		mX(c) = sqrt(ms/n);

      }
      else
      {
        //
        mX(c)=0;
        std::cout << "NaN column at " << c << "!" <<std::endl;
      }
    }

  }

  return mX;
}


ArrayType meanMat(Array2DType X, Array2DType P, int nanVal=0)
//take the column-wise 'weighted mean' of the matrix X, ignoring the zero elements.
{
  X = element_product(X, P);   //x = x.p(x)
  ArrayType mX;
  mX.SetSize(X.cols());
  ArrayType aCol;
  for (unsigned int c = 0; c<X.cols(); c++)
  {

    aCol = X.get_column(c);
    VariableType s;

    if (aCol.min_value() != nanVal)
    {
      s = aCol.sum();
      mX(c) = s/P.get_column(c).sum();
    }
    else
    {
      VariableType n = 0;
      for (unsigned int i =0; i<aCol.Size(); i++)
      {
        if (aCol(i) != nanVal)
        {
          n+= P.get_column(c).get(i);
          s+= aCol(i);
        }

      }
      if (n!=0)
      {
        mX(c) = s/n;
      }
      else
      {
        mX(c)=0;
        std::cout << "NaN column at " << c << "!" <<std::endl;
      }
    }

  }

  return mX;
}
VariableType Gamma(VariableType x, VariableType alpha, VariableType beta)
{
  VariableType gammaPdf;

  gammaPdf = 1/(vnl_gamma(alpha)*pow(beta,alpha))*pow(x,(alpha-1))*exp(-x/beta);

  return gammaPdf;
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

CurveType SmoothCurve(CurveType Curve)
{
  CurveType SmoothedCurve;
  int NumberOfPoints = Curve.rows();
  SmoothedCurve.set_size(NumberOfPoints,3);
  int window = 5;  //radius
  for (int j=0; j<NumberOfPoints; ++j)
  {
    CurvePointType sumPoints;
    sumPoints.set_size(3);
    sumPoints.fill(0);
    int el =0,low,up;
    if  ((j-window)<0)
    {
      low=0;
    }
    else
    {
      low =(j-window);
    }

    if  (NumberOfPoints<(j+window)) up=NumberOfPoints; else up =(j+window);
    for (int i=low; i< up; ++i)
    {
      sumPoints += Curve.get_row(i); el++;
    }
    SmoothedCurve.set_row(j, sumPoints/el);
  }
  //std::cout << SmoothedCurve <<std::endl;
  return SmoothedCurve;
}

CurveType SmoothAndResampleCurve(CurveType Curve, float ds)
{
	unsigned int numberOfLandmarks = Curve.rows();
	std::vector<CoordinateType> s = getArcLengthParameterization(Curve);

	typedef itk::ThinPlateSplineKernelTransform<double, 3 >  TransformType;
	TransformType::Pointer transform = TransformType::New();

	typedef TransformType::PointSetType       PointSetType;

	PointSetType::Pointer sourceLandmarks = PointSetType::New();
	PointSetType::Pointer targetLandmarks = PointSetType::New();

	transform->SetSourceLandmarks( sourceLandmarks );
	transform->SetTargetLandmarks( targetLandmarks );

	typedef PointSetType::PointsContainer 	PointsContainer;

	PointsContainer::Pointer sources = sourceLandmarks->GetPoints();
	PointsContainer::Pointer targets = targetLandmarks->GetPoints();

	transform->SetSourceLandmarks( sourceLandmarks );
	transform->SetTargetLandmarks( targetLandmarks );

	sources->Reserve( numberOfLandmarks );
	targets->Reserve( numberOfLandmarks );

	typedef PointSetType::PointType       PointType;
	PointType source;
	PointType target;

	for( unsigned int i = 0; i < numberOfLandmarks; i++ )
		{
		 target[0] = Curve.get(i,0);
		 target[1] = Curve.get(i,1);
		 target[2] = Curve.get(i,2);

		 //source[0] = static_cast<float>( i )/(numberOfLandmarks-1);
		 source[0] = s.at(i)/s.at(s.size()-1);
		 source[1] = 0;
		 source[2] = 0;
	     sources->InsertElement( i, source );
		 targets->InsertElement( i, target );
		}

  transform->SetStiffness(0);
  //A stiffness of zero results in the standard interpolating spline.
  //non-zero stiffness allows the spline to approximate rather than interpolate the landmarks.
  //Stiffness values are usually rather small, typically in the range of 0.001 to 0.1.
  transform -> UpdateParameters();

  transform->ComputeWMatrix();
  //resample:
  PointType temp_point, temp_outputpoint;
  CurveType smoothedCurve;
  unsigned int NumberOfPoints = static_cast<unsigned int> (s.at(s.size()-1)/ds);
  smoothedCurve.set_size(NumberOfPoints,3);

  for (unsigned int i = 0; i< NumberOfPoints; i++)   //check the size of vectors
  {
	  temp_point[0]= i*ds/s.at(s.size()-1); temp_point[1] =0; temp_point[2]=0;
	  //std::cout << temp_point << std::endl;
	  temp_outputpoint = transform->TransformPoint(temp_point);
	  //std::cout << temp_outputpoint << std::endl;
	  smoothedCurve.set_row(i, temp_outputpoint.GetVnlVector() );
  }

  return smoothedCurve;
}

std::vector<CoordinateType> getArcLengthParameterization(const CurveType curve)
{
	std::vector<CoordinateType> s;
	CurvePointType currentPoint, nextPoint;
	s.push_back(0); CoordinateType temp=0;
	for (unsigned int i=0; i< curve.rows()-1; i++)
	{
	    currentPoint = curve.get_row(i);
	    nextPoint = curve.get_row(i+1);

	    temp+= (nextPoint - currentPoint).two_norm();
	    s.push_back(temp);
	}

	return s;
}

VariableType diffCurve(CurveType MyCurve1, CurveType MyCurve2) //simplest implementation
{
  VariableType dist = 0;
  CurvePointType p1,p2;
  for (unsigned int l = 0; l< MyCurve2.rows(); ++l)
  {
    p2 = MyCurve2.get_row(l);
    p1 = MyCurve1.get_row(l);
    dist+=(p1 - p2).two_norm();
  }
  return dist/MyCurve2.rows();
}