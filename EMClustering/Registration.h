
#ifndef __Registration_h_included_
#define __Registration_h_included_


#include "Common.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAffineTransform.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCommand.h"

typedef itk::AffineTransform<CoordinateType,PointDimension>            TransformType;
typedef itk::RegularStepGradientDescentOptimizer                       OptimizerType;
typedef itk::MeanSquaresImageToImageMetric<ImageType,ImageType>        MetricType;
typedef itk::LinearInterpolateImageFunction<ImageType,CoordinateType>  InterpolatorType;


TransformType::Pointer doAffineRegistration(ImageType* , ImageType* , std::string);
MeshType::Pointer applyTransform(MeshType* , TransformType* , std::vector<unsigned long int> );

#endif // #ifndef
