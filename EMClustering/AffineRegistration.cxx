
#include <AffineRegistration.h>

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  RegistrationInterfaceCommand() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;
  void Execute(itk::Object * object, const itk::EventObject & event)
  {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }
    RegistrationPointer registration =
                        dynamic_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(
                       registration->GetOptimizer() );

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel()  << std::endl;
    std::cout << std::endl;

    if ( registration->GetCurrentLevel() == 0 )
      {
      optimizer->SetMaximumStepLength( 8.00 );
      optimizer->SetMinimumStepLength(  0.1 );
      }
    else
      {
      optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 4.0 );
      optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
      }
  }
  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }
};


class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};
  itk::ProcessObject::Pointer m_Registration;
  
public:
  typedef itk::RegularStepGradientDescentOptimizer  OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void SetRegistration( itk::ProcessObject *p)
    {
      m_Registration = p;
    }
  
  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer = 
        dynamic_cast< OptimizerPointer >( object );
      if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
        return;
        }
      
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetCurrentStepLength() << "   ";
      std::cout << optimizer->GetValue() << std::endl;
      if (m_Registration)
        {
        m_Registration->UpdateProgress( 
          static_cast<double>(optimizer->GetCurrentIteration()) /
          static_cast<double>(optimizer->GetNumberOfIterations()));
        }
    }
};


TransformType::Pointer doSlicerFastAffineRegistration(ImageType* fixedImage, ImageType* movingImage, std::string OutputDirectory)
{
  typedef itk::OrientImageFilter<ImageType,ImageType> FixedOrientFilterType;
  typedef itk::OrientImageFilter<ImageType,ImageType> MovingOrientFilterType;
  
  typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType>    MetricType;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateType>  InterpolatorType;
  //typedef itk::ImageRegistrationMethod<ImageType,ImageType>  RegistrationType;
  typedef OptimizerType::ScalesType OptimizerScalesType;
  typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleType;
  typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateType> ResampleInterpolatorType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef itk::ContinuousIndex<CoordinateType, 3> ContinuousIndexType;

	typedef itk::MultiResolutionImageRegistrationMethod< ImageType, ImageType    > RegistrationType;
	typedef itk::RecursiveMultiResolutionPyramidImageFilter<ImageType, ImageType  >  FixedImagePyramidType;
	typedef itk::RecursiveMultiResolutionPyramidImageFilter<ImageType, ImageType  >   MovingImagePyramidType;

	FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
	MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();


  FixedOrientFilterType::Pointer orientFixed = FixedOrientFilterType::New();//##
  //itk::PluginFilterWatcher watchOrientFixed(orientFixed,   "Orient Fixed Image",  CLPProcessInformation,  1.0/5.0, 0.0);
  orientFixed->UseImageDirectionOn();
  orientFixed->SetDesiredCoordinateOrientationToAxial();

  orientFixed->SetInput (fixedImage);

  orientFixed->Update();
  
  MovingOrientFilterType::Pointer orientMoving = MovingOrientFilterType::New();//##
  //itk::PluginFilterWatcher watchOrientMoving(orientMoving,  "Orient Moving Image", CLPProcessInformation,  1.0/5.0, 2.0/5.0);
  orientMoving->UseImageDirectionOn();
  orientMoving->SetDesiredCoordinateOrientationToAxial();
  
  orientMoving->SetInput (movingImage);
  
  orientMoving->Update();

  // Set up the optimizer
  //
  //
  int TranslationScale = 100;

    OptimizerType::Pointer      optimizer     = OptimizerType::New();
    optimizer->SetNumberOfIterations ( 1000 );
    optimizer->SetMinimumStepLength ( .0005 );
    optimizer->SetMaximumStepLength ( 10.0 );
    optimizer->SetMinimize(true);   

  TransformType::Pointer transform = TransformType::New();
  OptimizerScalesType scales( transform->GetNumberOfParameters() );
  scales.Fill ( 1.0 );
  for( unsigned j = 9; j < 12; j++ )
    {
    scales[j] = 1.0 / vnl_math_sqr(TranslationScale);
    }
  optimizer->SetScales( scales );

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  
  // Initialize the transform
  //
  TransformType::InputPointType centerFixed;
  ImageType::RegionType::SizeType sizeFixed = orientFixed->GetOutput()->GetLargestPossibleRegion().GetSize();
  // Find the center
  ContinuousIndexType indexFixed;
  for ( unsigned j = 0; j < 3; j++ )
    {
    indexFixed[j] = (sizeFixed[j]-1) / 2.0;
    }
  orientFixed->GetOutput()->TransformContinuousIndexToPhysicalPoint ( indexFixed, centerFixed );

  TransformType::InputPointType centerMoving;
  ImageType::RegionType::SizeType sizeMoving = orientMoving->GetOutput()->GetLargestPossibleRegion().GetSize();
  // Find the center
  ContinuousIndexType indexMoving;
  for ( unsigned j = 0; j < 3; j++ )
    {
    indexMoving[j] = (sizeMoving[j]-1) / 2.0;
    }
  orientMoving->GetOutput()->TransformContinuousIndexToPhysicalPoint ( indexMoving, centerMoving );

  transform->SetCenter( centerFixed );
  transform->Translate(centerMoving-centerFixed);
  std::cout << "Centering transform: "; //transform->Print( std::cout );
  
  // Set up the metric
  //
  int HistogramBins = 30;
  int SpatialSamples = 10000;
  MetricType::Pointer  metric   = MetricType::New();
  metric->SetNumberOfHistogramBins ( HistogramBins );
  metric->SetNumberOfSpatialSamples( SpatialSamples );
  metric->ReinitializeSeed(123);
  
  // Create the interpolator
  //
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  // Set up the registration
  //
  RegistrationType::Pointer registration = RegistrationType::New();
  registration->SetTransform ( transform );
  registration->SetInitialTransformParameters ( transform->GetParameters() );
  registration->SetMetric ( metric );
  registration->SetOptimizer ( optimizer );
  registration->SetInterpolator ( interpolator );
  registration->SetFixedImage ( orientFixed->GetOutput() );
  registration->SetMovingImage ( orientMoving->GetOutput() );

  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );

  registration->SetFixedImageRegion( orientFixed->GetOutput()->GetBufferedRegion() );
  registration->SetNumberOfLevels( 3 );

  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( itk::IterationEvent(), command );
  
  try
    {
    registration->Update();     
    } 
  catch( itk::ExceptionObject & err )
    {
    std::cout << err << std::endl;
    std::cerr << err << std::endl;
   // return  EXIT_FAILURE ;
    } 
  catch ( ... )
    {
   // return  EXIT_FAILURE ;
    }

  transform->SetParameters ( registration->GetLastTransformParameters() );

  // Resample to the original coordinate frame (not the reoriented
  // axial coordinate frame) of the fixed image
  //
    ResampleType::Pointer resample = ResampleType::New();
    ResampleInterpolatorType::Pointer Interpolator = ResampleInterpolatorType::New();
    
    resample->SetInput ( movingImage );
    resample->SetTransform ( transform );
    resample->SetInterpolator ( Interpolator );

    // Set the output sampling based on the fixed image.
    // ResampleImageFilter needs an image of the same type as the
    // moving image.
    ImageType::Pointer fixedInformation = ImageType::New();
    fixedInformation->CopyInformation( fixedImage );
    resample->SetOutputParametersFromImage ( fixedInformation );

    resample->Update();
    std::cout << "Writing the transformed FA volume ..." << std::endl;
    std::string filename = OutputDirectory + "/TransformedFA.nhdr";
    
    WriterType::Pointer resampledWriter = WriterType::New();
    resampledWriter->SetFileName ( filename );
    resampledWriter->SetInput ( resample->GetOutput() );
    try
      {
      resampledWriter->Write();
      }
    catch( itk::ExceptionObject & err )
      { 
      std::cerr << err << std::endl;
      std::cerr << err << std::endl;
      //return EXIT_FAILURE;
      }
  return transform;
}
MeshType::Pointer applyTransform(MeshType* atlasCenters, TransformType* transform, std::vector<unsigned long int> CellIDs)
{
	//inverse affine transform should be applied to the points:
	//transform->Print(std::cout);
    //TransformType::ParametersType invTransParams;
	TransformType::Pointer invTransform = TransformType::New();
	invTransform->SetCenter(transform->GetCenter());
	transform->GetInverse(invTransform);
	//invTransform->Print(std::cout);

	MeshType::Pointer transformedCenters=MeshType::New();
	CellAutoPointer aCell, MyCell;
	long int myid =0;

	for (unsigned int c=0; c<CellIDs.size(); ++c)
	{
		MyCell.TakeOwnership( new PolylineType );
		MeshType::CellIdentifier CellID = CellIDs.at(c);
		atlasCenters->GetCell(CellID, aCell);
		CellDataType cellvalue;
		atlasCenters->GetCellData(CellID, &cellvalue);
		PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
		MeshType::PointType point, tpoint;
		MeshType::PixelType pointvalue;

		for (unsigned int j=0; j < aCell->GetNumberOfPoints(); j++)
		{
			atlasCenters->GetPoint(*pit, &point);
			atlasCenters->GetPointData(*pit, &pointvalue);
			tpoint = invTransform->TransformPoint(point);
			transformedCenters->SetPoint(myid, tpoint );
			transformedCenters->SetPointData(myid, pointvalue );
			MyCell->SetPointId(j,myid);
			pit++; myid++;
		}
		transformedCenters->SetCell(c, MyCell);
		transformedCenters->SetCellData(c, cellvalue);
	}

	return transformedCenters;
}

CenterType applyTransform(CenterType meshes, TransformType* transform)
{

	CenterType transformed_meshes;
	std::vector<unsigned long int> CellIDs;
	CellIDs.push_back(0);
	for (unsigned int k=0; k<meshes.size(); k++)
	{
		MeshType::Pointer t = applyTransform(meshes.at(k), transform, CellIDs);
		transformed_meshes.push_back(t);
	}

	return transformed_meshes;
}

  
