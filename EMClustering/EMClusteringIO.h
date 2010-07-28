
#ifndef __EMClusteringIO_h_included_
#define __EMClusteringIO_h_included_

#include "Common.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
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


void WriteVTKfile(MeshType*, std::string, CopyFieldType);

MeshType::Pointer ReadVTKfile(std::string);

MeshType::Pointer ReadVTKfiles(std::vector<std::string>);

ImageType::Pointer ReadImageVolume(std::string);

void WriteImageVolume(ImageType* , std::string);

void WriteCSVfile(std::string, const Array2DType &);

void writeMCSVfile(std::string, const ArrayType, const ArrayType, const std::vector<std::string> & );

#endif // #ifndef
