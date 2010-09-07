
#include <myMaths.h>

ArrayType meanMat(const Array2DType &X, int nanVal=0)
//take the column-wise mean of the matrix X, ignoring the zero elements.
{
  std::vector<int> nanValCols;
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
        nanValCols.push_back(c);
      }
     }

  }
  std::cout << "NaN colums: "<< nanValCols.size()<<std::endl;
  return mX;
}

ArrayType stdMat(const Array2DType &X, int nanVal=0)
//take the column-wise std of the matrix X, ignoring the nonVal elements.
{
  std::vector<int> nanValCols;
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
        nanValCols.push_back(c);
        }
      }

  }
  std::cout << "NaN colums: "<< nanValCols.size()<<std::endl;
  return mX;
}

ArrayType meanMat(Array2DType X, Array2DType P, int nanVal=0)
//take the column-wise 'weighted mean' of the matrix X, ignoring the zero elements.
{
  std::vector<int> nanValCols;
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
        nanValCols.push_back(c);
      }
    }

  }
  std::cout << "NaN colums: "<< nanValCols.size()<<std::endl;
  return mX;
}
VariableType Gamma(VariableType x, VariableType alpha, VariableType beta)
{
  VariableType gammaPdf;

  gammaPdf = 1/(vnl_gamma(alpha)*pow(beta,alpha))*pow(x,(alpha-1))*exp(-x/beta);

  return gammaPdf;
}



CurveType SmoothCurve(CurveType Curve, VariableType samplesDistance)
{
  CurveType SmoothedCurve;
  int NumberOfPoints = Curve.rows();
  SmoothedCurve.set_size(NumberOfPoints,3);
  VariableType windowLength = 6;  // in mm -- parameter
  int  window = ceil(windowLength/2/samplesDistance);

  SmoothedCurve.set_row (0, Curve.get_row(0));
  for (int j=1; j<NumberOfPoints-1; ++j)
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
  SmoothedCurve.set_row (NumberOfPoints-1, Curve.get_row(NumberOfPoints-1));
  //std::cout << SmoothedCurve <<std::endl;
  return SmoothedCurve;
}

CurveType SmoothAndResampleCurve(CurveType Curve, VariableType ds)
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

  transform->SetStiffness(0 );
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

  for (unsigned int i = 0; i< NumberOfPoints; i++)
  {
	  temp_point[0]= static_cast<float> (i)/(NumberOfPoints-1); temp_point[1] =0; temp_point[2]=0;
	  temp_outputpoint = transform->TransformPoint(temp_point);
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
  unsigned int numberOfRows = MyCurve2.rows();
  if (MyCurve2.rows()>MyCurve1.rows())
  {
	  numberOfRows = MyCurve1.rows();
  }
  for (unsigned int l = 0; l< numberOfRows; ++l)
  {
    p2 = MyCurve2.get_row(l);
    p1 = MyCurve1.get_row(l);
    dist+=(p1 - p2).two_norm();
  }
  return dist/MyCurve2.rows();
}
