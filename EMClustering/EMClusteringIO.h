#include "itksys/SystemTools.hxx"
#include "itksys/Glob.hxx"
#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkPolylineCell.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkArray2D.h>
#include <itkDiffusionTensor3D.h>


#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkCellArray.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkCellData.h"
#include <vtkDataSetAttributes.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkLongArray.h>

const unsigned int PointDimension = 3;
const unsigned int MaxTopologicalDimension = 1;
typedef double     CoordinateType;
typedef float      VariableType;
typedef itk::Vector<CoordinateType, PointDimension> VectorType;
typedef itk::Array<VariableType>   ArrayType;
typedef itk::Array2D<VariableType> Array2DType;
typedef std::vector<Array2DType>   Array3DType;

typedef struct{double                      FA;
itk::FixedArray<double, PointDimension >   EigenValues;
itk::FixedArray<double, 9 >                Tensor;
itk::Point<CoordinateType, PointDimension> AtlasPosition;
itk::Array<long int>                       Correspondence;
}PixelType;

typedef struct{std::string   CaseName;
int                          ClusterLabel;
ArrayType                    membershipProbability;
ArrayType                    atlasPriors;
}CellDataType;

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




void WriteVTKfile(MeshType*, std::string, CopyFieldType);

MeshType::Pointer ReadVTKfile(std::string);

MeshType::Pointer ReadVTKfiles(std::vector<std::string>);

void WriteCSVfile(std::string, const Array2DType &);

void writeMCSVfile(std::string, const ArrayType, const ArrayType, const std::vector<std::string> & );
