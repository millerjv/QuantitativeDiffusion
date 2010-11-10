
#ifndef __myMaths_h_included_
#define __myMaths_h_included_

#include "Common.h"
#include "vnl/vnl_gamma.h"
#include <itkThinPlateSplineKernelTransform.h>

typedef itk::Array2D<CoordinateType>                        CurveType;
typedef itk::Array<CoordinateType>                          CurvePointType;


ArrayType meanMat(Array2DType, int);
ArrayType stdMat(Array2DType, int );
ArrayType meanMat(Array2DType, Array2DType, int);
VariableType Gamma(VariableType, VariableType, VariableType);
CurveType SmoothCurve(CurveType, VariableType);
VariableType diffCurve(CurveType, CurveType);
CurveType SmoothAndResampleCurve(CurveType, VariableType);
std::vector<CoordinateType> getArcLengthParameterization(CurveType);

#endif // #ifndef
