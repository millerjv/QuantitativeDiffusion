
#ifndef __AffineRegistration_h_included_
#define __AffineRegistration_h_included_


#include "Common.h"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "itkOrientedImage.h"
#include "itkOrientImageFilter.h"

#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkBinomialBlurImageFilter.h"
#include "itkCommand.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


typedef itk::AffineTransform<CoordinateType>            TransformType;

// The function is written similarly to Slicer fast affine registration module:
TransformType::Pointer doSlicerFastAffineRegistration(ImageType* , ImageType*, std::string);
MeshType::Pointer applyTransform(MeshType* , TransformType* , std::vector<unsigned long int> );
CenterType applyTransform(CenterType, TransformType*);

#endif // #ifndef
