
#ifndef __Quantification_h_included_
#define __Quantification_h_included_

#include "Common.h"

void ComputeScalarMeasures(MeshType*);
Array3DType BuildFeatureMatrix(const MeshType*, const MeshType*, int);
Array3DType BuildFeatureMatrix(const MeshType*, const QuadEdgeMeshType*, int);

#endif // #ifndef
