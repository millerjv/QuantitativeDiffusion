
#ifndef __Common_h_included_
#define __Common_h_included_

#include "itksys/SystemTools.hxx"
#include "itksys/Glob.hxx"
#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkPolylineCell.h"
#include "itkDefaultStaticMeshTraits.h"
#include <itkArray2D.h>
#include <itkDiffusionTensor3D.h>
#include "itkOrientedImage.h"
#include "itkImage.h"

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"

const unsigned int PointDimension = 3;
const unsigned int MaxTopologicalDimension = 1;
typedef double     CoordinateType;
typedef float      VariableType;
typedef itk::Vector<CoordinateType, PointDimension> VectorType;
typedef itk::Array<VariableType>   ArrayType;
typedef itk::Array2D<VariableType> Array2DType;
typedef std::vector<Array2DType>   Array3DType;
typedef itk::OrientedImage<CoordinateType, PointDimension>  ImageType;

typedef struct{double                      FA;
itk::FixedArray<double, PointDimension >   EigenValues;
itk::FixedArray<double, 9 >                Tensor;
itk::Point<CoordinateType, PointDimension> AtlasPosition;
itk::Array<long int>                       Correspondence;
}PixelType;

typedef struct{
itk::FixedArray<double, PointDimension >   Orientation;
}QEPixelType;


typedef struct{std::string   CaseName;
int                          ClusterLabel;
ArrayType                    membershipProbability;
ArrayType                    atlasPriors;
}CellDataType;

typedef struct{
int                          ClusterLabel;
}QECellDataType;

typedef struct{bool   FA;
bool                  EigenValues;
bool                  Tensor;
bool                  ClusterLabel;
bool                  CaseName;
bool                  Correspondences;
}CopyFieldType;


typedef double InterpolationWeightType;
typedef itk::DefaultStaticMeshTraits<
PixelType, PointDimension, MaxTopologicalDimension,
CoordinateType, InterpolationWeightType, CellDataType >     MeshTraits;
typedef itk::Mesh< PixelType, PointDimension, MeshTraits >  MeshType;
typedef MeshType::CellType                                  CellType;
typedef itk::PolylineCell< CellType >                       PolylineType;
typedef CellType::CellAutoPointer                           CellAutoPointer;
typedef itk::DiffusionTensor3D< CoordinateType >            TensorPixelType;
typedef TensorPixelType::RealValueType                      RealValueType;
typedef TensorPixelType::EigenValuesArrayType               EigenValuesArrayType;

typedef itk::QuadEdgeMeshTraits<
QEPixelType, PointDimension, QEPixelType, QECellDataType,
CoordinateType, InterpolationWeightType >     QuadEdgeMeshTraits;

typedef itk::QuadEdgeMesh< QEPixelType, PointDimension, QuadEdgeMeshTraits>  QuadEdgeMeshType;
typedef std::vector<QuadEdgeMeshType::Pointer>              CenterType;

#endif // #ifndef
