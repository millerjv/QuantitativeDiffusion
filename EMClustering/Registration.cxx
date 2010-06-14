
#include "EMClusteringIO.h"
#include "Registration.h"

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
                      dynamic_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition();
    }
};

TransformType::Pointer doAffineRegistration(ImageType* caseFAVolume, ImageType* atlasFAVolume)
{
  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  TransformType::Pointer  transform = TransformType::New();
  registration->SetTransform( transform );

  registration->SetFixedImage(caseFAVolume );
  registration->SetMovingImage(atlasFAVolume);

  registration->SetFixedImageRegion( caseFAVolume->GetBufferedRegion() );

  typedef itk::CenteredTransformInitializer< TransformType, ImageType, ImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage( caseFAVolume );
  initializer->SetMovingImage( atlasFAVolume );
  initializer->MomentsOn();
  initializer->InitializeTransform();
  registration->SetInitialTransformParameters( transform->GetParameters() );
  double translationScale = 1.0 / 1000.0;
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  1.0;
  optimizerScales[5] =  1.0;
  optimizerScales[6] =  1.0;
  optimizerScales[7] =  1.0;
  optimizerScales[8] =  1.0;
  optimizerScales[9]  =  translationScale;
  optimizerScales[10] =  translationScale;
  optimizerScales[11] =  translationScale;

  optimizer->SetScales( optimizerScales );
  unsigned int maxNumberOfIterations = 10;
  optimizer->SetMaximumStepLength( 0.1 );
  optimizer->SetMinimumStepLength( 0.001 );
  optimizer->SetNumberOfIterations( maxNumberOfIterations );
  optimizer->MinimizeOn();

  //CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  //optimizer->AddObserver( itk::IterationEvent(), observer );

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    //return EXIT_FAILURE;
    }
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetCenter( transform->GetCenter() );
  finalTransform->SetParameters( finalParameters );

  return finalTransform;
}

MeshType::Pointer applyTransform(MeshType* atlasCenters, TransformType* transform, std::vector<unsigned long int> CellIDs)
{

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
    TransformType::ParametersType invTransParams;
    TransformType::Pointer invTransform = TransformType::New();
    invTransParams = transform->GetInverseTransform()->GetParameters();
    invTransform->SetCenter( transform->GetCenter() );
    invTransform->SetParameters( invTransParams );
    for (unsigned int j=0; j < aCell->GetNumberOfPoints(); j++)
    {
      atlasCenters->GetPoint(*pit, &point);
      atlasCenters->GetPointData(*pit, &pointvalue);
      //tpoint = transform->GetInverseTransform()->TransformPoint(point);
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

