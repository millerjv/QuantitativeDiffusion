#ifndef __EMCore_h_included_
#define __EMCore_h_included_

#include "Common.h"

void  setPriorInfo(Array2DType &, MeshType* );
void SetInitialValue(const Array2DType &, ArrayType &);
ArrayType AdjustThreshold(VariableType, ArrayType, ArrayType);
Array2DType ComputePosterior(const Array2DType &, const Array2DType &);
Array2DType ComputeLikelihood(const Array2DType &, ArrayType, ArrayType);
void UpdateModelParameters(const Array2DType &, const Array2DType &, ArrayType& , ArrayType& , Array2DType& , bool);
void AssignClusterLabels(MeshType*, const Array2DType &);


#endif // #ifndef
