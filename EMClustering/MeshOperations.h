
#ifndef __MeshOperations_h_included_
#define __MeshOperations_h_included_

#include <Common.h>
#include "myMaths.h"
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkBoundingBox.h>
#include <itkLinearInterpolateImageFunction.h>

Array2DType ComputeDissimilarity(MeshType*, MeshType*, ImageType*, VariableType, VariableType, VariableType);
MeshType::Pointer UpdateCenters(MeshType*, MeshType*, const Array2DType &, VariableType);
MeshType::Pointer RefineData(const MeshType*, Array2DType &, Array2DType &, Array2DType &, ArrayType, bool);
MeshType::Pointer SmoothMesh(MeshType*, VariableType, bool);
ArrayType diffMeshes(const MeshType*, const MeshType*);
MeshType::Pointer getTrajectories(MeshType*, std::vector<unsigned long int>);
MeshType::Pointer getCluster(MeshType*,int);
MeshType::Pointer getCluster(MeshType*, std::string);
Array2DType getClusterPosterior(Array2DType, MeshType*,int);
std::vector<std::string> getClusterSubjectNames(MeshType*);
VariableType getSampleSpacing(const MeshType*);
ImageType::Pointer getSubSpace(const MeshType*, VariableType);
void  AddPointScalarToACell(MeshType*, MeshType::CellIdentifier, ArrayType);


#endif // #ifndef
