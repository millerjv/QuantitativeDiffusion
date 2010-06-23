
#include "vnl/vnl_gamma.h"

typedef itk::Array2D<CoordinateType>                        CurveType;
typedef itk::Array<CoordinateType>                          CurvePointType;


ArrayType meanMat(const Array2DType &, int );
ArrayType stdMat(const Array2DType &, int );
ArrayType meanMat(Array2DType, Array2DType, int);
VariableType Gamma(VariableType, VariableType, VariableType);
Array2DType ComputePosterior(const Array2DType &, const Array2DType &);
Array2DType ComputeLikelihood(const Array2DType &, ArrayType, ArrayType);
CurveType SmoothCurve(CurveType);
VariableType diffCurve(CurveType, CurveType);
CurveType SmoothAndResampleCurve(CurveType, float);
std::vector<CoordinateType> getArcLengthParameterization(CurveType);
