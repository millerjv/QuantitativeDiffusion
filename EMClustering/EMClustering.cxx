#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
//// SLICER VERSION
#include "itksys/SystemTools.hxx"
#include "itksys/Glob.hxx"
#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkPolylineCell.h"
#include "itkDefaultStaticMeshTraits.h"
#include <itkDiffusionTensor3D.h> 
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkOrientedImage.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include <itkMatrix.h>
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkArray2D.h>
#include <itkBoundingBox.h>
#include <itkExtractImageFilter.h>
#include "vnl/vnl_gamma.h"
#include "EMClusteringCLP.h"

#include "vtkActor.h"
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

#include "itkAffineTransform.h"
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
typedef itk::OrientedImage<CoordinateType,PointDimension >  ImageType;
typedef itk::Array2D<CoordinateType>                        CurveType;
typedef itk::Array<CoordinateType>                          CurvePointType;
typedef itk::AffineTransform<CoordinateType,PointDimension>            TransformType;
typedef itk::RegularStepGradientDescentOptimizer                       OptimizerType;
typedef itk::MeanSquaresImageToImageMetric<ImageType,ImageType>        MetricType;
typedef itk::LinearInterpolateImageFunction<ImageType,CoordinateType>  InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType,ImageType >             RegistrationType;


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

      /* Print the angle for the trace plot
      vnl_matrix<double> p(2, 2);
      p[0][0] = (double) optimizer->GetCurrentPosition()[0];
      p[0][1] = (double) optimizer->GetCurrentPosition()[1];
      p[1][0] = (double) optimizer->GetCurrentPosition()[2];
      p[1][1] = (double) optimizer->GetCurrentPosition()[3];
      vnl_svd<double> svd(p);
      vnl_matrix<double> r(2, 2);
      r = svd.U() * vnl_transpose(svd.V());
      double angle = asin(r[1][0]);
      std::cout << " AffineAngle: " << angle * 45.0 / atan(1.0) << std::endl;
      */
    }
};

void CopyItkMesh2VtkPolyData(MeshType* mesh, vtkPolyData* polydata, CopyFieldType copyField)

// Convert the itk mesh to vtk polydata:
{
  unsigned int numPoints = mesh->GetNumberOfPoints();
  unsigned int numCells  = mesh->GetNumberOfCells();

  vtkPoints* vpoints = vtkPoints::New();
  vpoints->SetNumberOfPoints(numPoints);

  vtkDoubleArray* tensors = vtkDoubleArray::New();
  tensors->SetNumberOfComponents(9);
  tensors->SetNumberOfTuples(numPoints);

  vtkDoubleArray* scalars = vtkDoubleArray::New();
  scalars->SetNumberOfTuples(numPoints);
  scalars->SetName("FA");

  vtkUnsignedLongArray* clusterScalars = vtkUnsignedLongArray::New();
  clusterScalars->SetNumberOfTuples(numCells);
  clusterScalars->SetName("ClusterId");


  vtkDoubleArray* clusterMembershipProbs = vtkDoubleArray::New();
  clusterMembershipProbs->SetNumberOfTuples(numCells);
  clusterMembershipProbs->SetName("membershipProbabilities");

  vtkStringArray* subjectName = vtkStringArray::New();
  subjectName->SetName("CaseName");
  subjectName->SetNumberOfTuples(numCells);

  itk::FixedArray<double, 9 >  MyTensor;
  double MyFA;
  MeshType::PixelType pointvalue;
  MeshType::PointsContainer::Pointer points = mesh->GetPoints();
  for(MeshType::PointsContainer::Iterator i = points->Begin(); i !=
    points->End(); ++i)
  {
    int idx = i->Index();
    MeshType::PointType ip = i->Value();
    //take care of orientation difference between itk and vtk:
    ip[0] = -ip[0];
    ip[1] = -ip[1];

    vpoints->SetPoint(idx, ip[0], ip[1], ip[2]);
    mesh->GetPointData(idx, &pointvalue);
    MyTensor = pointvalue.Tensor;
    MyFA = pointvalue.FA;
    scalars->InsertTuple1(idx, MyFA);
    tensors->SetTuple9(idx,MyTensor[0],MyTensor[1],MyTensor[2],MyTensor[3],MyTensor[4],MyTensor[5],MyTensor[6],MyTensor[7],MyTensor[8]);
  }

  polydata->SetPoints(vpoints);

  //COPY POINT DATA

  if (copyField.Tensor)
  {
    polydata->GetPointData()->SetTensors(tensors);
  }

  if (copyField.FA)
  {
    polydata->GetPointData()->SetScalars(scalars);
  }


  vtkCellArray *polylines = vtkCellArray::New();
  CellAutoPointer acell;
  CellDataType cellvalue;
  for (unsigned int i=0; i < mesh->GetNumberOfCells(); ++i)
  {
    mesh->GetCell(i, acell);
    polylines->InsertNextCell(acell->GetNumberOfPoints());
    PolylineType::PointIdIterator pit = acell->PointIdsBegin();
    for (unsigned int j=0; j < acell->GetNumberOfPoints(); ++j)
    {
      polylines->InsertCellPoint(*pit);
      ++pit;
    }

    mesh->GetCellData(i, &cellvalue);
    clusterScalars->SetValue(i,cellvalue.ClusterLabel);
    std::string filename = itksys::SystemTools::GetFilenameName(cellvalue.CaseName);
    std::string subfilename = filename.substr(0,13);
    //std::cout << subfilename.c_str()<< std::endl;
    subjectName->SetValue(i, subfilename.c_str()); //cellvalue.CaseName
    clusterMembershipProbs->SetNumberOfComponents(cellvalue.membershipProbability.Size());
    for (unsigned int p=0; p< cellvalue.membershipProbability.Size(); p++)
    {
      clusterMembershipProbs->InsertComponent(i,p,cellvalue.membershipProbability(p));
    }
  }

  //COPY CELL DATA
  if (copyField.ClusterLabel)
  {
    polydata->GetCellData()->AddArray(clusterScalars);
    polydata->GetCellData()->AddArray(clusterMembershipProbs);
  }

  if (copyField.CaseName)
  {
    polydata->GetCellData()->AddArray(subjectName);
  }
  polydata->SetLines( polylines );

  polylines->Delete();
  vpoints->Delete();
  tensors->Delete();
  scalars->Delete();
  clusterScalars->Delete();
  clusterMembershipProbs->Delete();
  subjectName->Delete();
}

vtkPolyData* itk2vtkPolydata(MeshType* mesh, CopyFieldType copyField)
{
  vtkPolyData* polydata = vtkPolyData::New();
  CopyItkMesh2VtkPolyData(mesh, polydata, copyField);
  return polydata;
}

void CopyVtkPolyData2ItkMesh(vtkPolyData* polydata, MeshType* mesh)
// Convert the vtk polydata to itk mesh:
{
  vtkPoints* vpoints = polydata->GetPoints();
  int numPoints = polydata->GetNumberOfPoints();
  vtkCellArray *polylines = polydata->GetLines();

  vtkDataArray *tensors = polydata->GetPointData()->GetTensors();
  vtkDataArray *clusterScalars = polydata->GetCellData()->GetScalars("ClusterId");

  MeshType::PixelType pointvalue;
  for(int i=0; i<numPoints; ++i)
  {
    // take care of the orientation difference between itk and vtk
    MeshType::PointType vpoint = vpoints->GetPoint(i);
    vpoint[0]= - vpoint[0];
    vpoint[1]= - vpoint[1];

    mesh->SetPoint(i, vpoint);
    if (tensors)
    {
      pointvalue.Tensor = tensors->GetTuple9(i);
      mesh->SetPointData(i, pointvalue);
    }
  }

  CellAutoPointer acell;
  CellDataType cellvalue;
  vtkIdType npts;
  vtkIdType *pts;
  polylines->InitTraversal();
  for (int j=0; j < polydata->GetNumberOfLines(); ++j)
  {
    acell.TakeOwnership( new PolylineType );
    polylines->GetNextCell(npts, pts);
    //acell->SetPointIds((unsigned long*)pts, (unsigned long *)&(pts[npts-1]));
    for (int jj=0; jj < npts; ++jj)
    {
      acell->SetPointId(jj, (CellType::PointIdentifier) pts[jj]);
    }
    mesh->SetCell(j, acell);
    if (clusterScalars)
    {
      cellvalue.ClusterLabel = clusterScalars->GetTuple1(j);
      mesh->SetCellData(j, cellvalue);

    }
  }

}
MeshType::Pointer vtk2itkMesh(vtkPolyData* polydata)
{
  MeshType::Pointer mesh = MeshType::New();
  CopyVtkPolyData2ItkMesh(polydata, mesh);
  return mesh;
}

void WriteVTKfile(MeshType* mesh, std::string filename, CopyFieldType copyField)

{
  vtkPolyData* polydata;
  polydata = itk2vtkPolydata(mesh, copyField);

  vtkXMLPolyDataWriter *MyPolyDataWriter = vtkXMLPolyDataWriter::New();
  MyPolyDataWriter->SetFileName( filename.c_str() );
  MyPolyDataWriter->SetInput(polydata);
  std::cout<< "Writing out "<< filename.c_str() <<"..."<<std::endl;
  MyPolyDataWriter->Update();
  MyPolyDataWriter->Delete();
  polydata->Delete();
}


void addMesh(MeshType* popMesh, MeshType* caseMesh, const std::string caseName)
{
  //add caseMesh to popMesh
	unsigned int pitNew = popMesh->GetNumberOfPoints();
	unsigned int citNew = popMesh->GetNumberOfCells();

	for(unsigned int i=0; i<caseMesh->GetNumberOfCells(); i++)
    {
	    CellAutoPointer atrajectory, copiedtrajectory;
	    copiedtrajectory.TakeOwnership( new PolylineType );

        caseMesh->GetCell(i, atrajectory);
        PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
        MeshType::PointType point;
        MeshType::PixelType pointvalue;
        CellDataType cellvalue;

        for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
        {
          caseMesh->GetPoint(*pit, &point);
          caseMesh->GetPointData( *pit, &pointvalue );
          popMesh->SetPoint(pitNew, point);
          popMesh->SetPointData( pitNew, pointvalue );
          copiedtrajectory->SetPointId(j,pitNew);
          ++pit; ++pitNew;
        }

        caseMesh->GetCellData(i, &cellvalue);
        popMesh->SetCell(citNew, copiedtrajectory);
        cellvalue.CaseName = caseName;
        popMesh->SetCellData(citNew, cellvalue);
        citNew++;
    }
}

MeshType::Pointer ReadVTKfile(std::string filename)
{
  MeshType::Pointer     mesh;
  std::string extension = itksys::SystemTools::GetFilenameExtension(filename);
  if (extension.compare(".VTP")==0 || extension.compare(".vtp")==0)
  {
	  vtkXMLPolyDataReader *MyPolyDataReader = vtkXMLPolyDataReader::New();
	  MyPolyDataReader->SetFileName( filename.c_str() );
	  std::cout<< "Reading "<<filename.c_str()<< "..."<<std::endl;
	  MyPolyDataReader->Update();
	  vtkPolyData* rpolydata = MyPolyDataReader->GetOutput();
	  //rpolydata->Print(std::cout);
      mesh = vtk2itkMesh(rpolydata);
      MyPolyDataReader->Delete();
  }
  else if (extension.compare(".VTK")==0 || extension.compare(".vtk")==0)
    {
  	  vtkPolyDataReader *MyPolyDataReader = vtkPolyDataReader::New();
  	  MyPolyDataReader->SetFileName( filename.c_str() );
  	  std::cout<< "Reading "<<filename.c_str()<< "..."<<std::endl;
  	  MyPolyDataReader->Update();
  	  vtkPolyData* rpolydata = MyPolyDataReader->GetOutput();
  	  mesh = vtk2itkMesh(rpolydata);
      MyPolyDataReader->Delete();
    }
  else
  {
	std::cerr<< extension << " is not a valid extension!" <<std::endl;
  }
  return mesh;

}

MeshType::Pointer ReadVTKfiles(std::vector<std::string> allfilenames)
{
	// create a new mesh
	MeshType::Pointer popMesh = MeshType::New();
	// for loop over the number of files
	for (unsigned int v=0; v< allfilenames.size(); v++)
	{
     // read each vtk file
	std::string caseName = allfilenames.at(v);
	MeshType::Pointer caseMesh =  ReadVTKfile(caseName);
	// add the existing mesh to the big mesh
	addMesh(popMesh, caseMesh, caseName);
	}

	return popMesh;
}

void WriteCSVfile(std::string fileName, const Array2DType &mat)
{
  ofstream myfile;
  std::cout<< "Writing out "<< fileName.c_str() << "..." << std::endl;
  myfile.open (fileName.c_str());
  for (unsigned long int r=0; r<mat.rows(); r++)
  {
    for (unsigned int c=0; c<mat.cols(); c++)
    {
      myfile << mat(r,c);
      if (c <mat.cols()-1)
      {
        myfile << ",";
      }
    }
    if (r<mat.rows()-1)
    {
      myfile << std::endl;
    }
  }
  myfile.close();
}

void writeMCSVfile(std::string fileName, const ArrayType &y, const ArrayType &yerr, const std::vector<std::string> &labels )
{
  ofstream myfile;
  std::cout<< "Writing out "<< fileName.c_str() << "..." << std::endl;
  myfile.open (fileName.c_str());
  unsigned int n = y.size();
  myfile << labels.at(0) << "," << labels.at(1) << "," << labels.at(2) << std::endl;
  if (n>1)
  {
   for (unsigned long int r=0; r<n; r++)
   {
     myfile << r << "," << y(r) << "," << yerr(r) <<std::endl;
   }
  }
  myfile.close();
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Array2DType ComputeDissimilarity(MeshType* mesh, MeshType* mesh_centers, ImageType* space)
{
  unsigned long int NumberOfTrajectories=mesh->GetNumberOfCells();
  unsigned int NumberOfClusters=mesh_centers->GetNumberOfCells();
  Array2DType DissimilarityMeasure;
  DissimilarityMeasure.SetSize(NumberOfTrajectories,NumberOfClusters);
  VariableType LargeDist = itk::NumericTraits<VariableType>::max();
  DissimilarityMeasure.fill(LargeDist);


  for (unsigned int ClusterIdx = 0; ClusterIdx<NumberOfClusters; ++ClusterIdx)
  {
    // Build a distance map for each cluster center
    space->FillBuffer(0);

    // Mark up the voxels of the centerlines on the image
    CellAutoPointer Centerline;
    mesh_centers->GetCell(ClusterIdx, Centerline);
    PolylineType::PointIdIterator pit = Centerline->PointIdsBegin();
    ImageType::IndexType ind;
    MeshType::PointType point;
    std::vector<long int> MyLabels;
    MyLabels.clear();
    bool outOfSpace = 1;
    for (unsigned int j=0; j < Centerline->GetNumberOfPoints(); ++j)
    {
      mesh_centers->GetPoint(*pit, &point);
      if (space->TransformPhysicalPointToIndex(point, ind))
      {
        outOfSpace = 0;
        space->SetPixel(ind,space->ComputeOffset(ind));
        //This part tries to populate the unique labels in MyLabels along each center
        if (MyLabels.size()==0)
        {
          MyLabels.push_back(space->ComputeOffset(ind));
        }
        else
        {
          if (!(MyLabels.at(MyLabels.size()-1)==space->ComputeOffset(ind)))
          {
            MyLabels.push_back(space->ComputeOffset(ind));
          }
        }
        //-------------
      }
      else
      {
        std::cout<<"Point "<< point<<" on the center "<< ClusterIdx <<" is out of the space"<<std::endl;
      }
      ++pit;
    }


    if (!outOfSpace)
    {
      std::cout<< " Generating the Distance Map for Cluster "<< ClusterIdx+1 <<" ..."<< std::endl;

      // Apply the daneilsson filter
      typedef itk::DanielssonDistanceMapImageFilter< ImageType, ImageType > DistanceMapFilterType;
      DistanceMapFilterType::Pointer DMFilter = DistanceMapFilterType::New();
      DMFilter->SetInput(space);
      DMFilter->InputIsBinaryOff();
      DMFilter->Update();

      ImageType::Pointer DistanceMap=DMFilter->GetOutput();
      ImageType::Pointer VoronoiMap=DMFilter->GetVoronoiMap();

      //write out the images:
      /*typedef itk::ImageFileWriter< ImageType > WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName("myspace.nhdr");
      writer->SetInput(VoronoiMap);
      std::cout<< "Writing the voronoi map" << std::endl;
      writer->Update(); */

      //Create the interpolator for the Distance Map
      itk::LinearInterpolateImageFunction<ImageType, double>::Pointer DMInterpolator =
        itk::LinearInterpolateImageFunction<ImageType, double>::New();
      DMInterpolator->SetInputImage(DistanceMap);

      //Find the dissimilarity measure for each trajectory
      for (unsigned int t=0; t<NumberOfTrajectories; ++t)
      {
        CellAutoPointer atrajectory;
        mesh->GetCell(t, atrajectory);
        PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
        ImageType::IndexType ind;
        MeshType::PointType point;
        MeshType::PixelType pointvalue;
        VariableType sumdist = 0;
        std::vector<long int> MyLabelsOnTrajectory;
        MyLabelsOnTrajectory.clear();

        for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
        {
          mesh->GetPoint(*pit, &point);
          if (space->TransformPhysicalPointToIndex(point, ind))
          {
            itk::ContinuousIndex<double, 3> cidx;
            if (space->TransformPhysicalPointToContinuousIndex( point, cidx ))
            {
              sumdist+=DMInterpolator->Evaluate(point);}
            mesh->GetPointData( *pit, &pointvalue );
            pointvalue.Correspondence.set_size(NumberOfClusters);
            //1st output -- the correspondence info is going to be used in updating
            //the centers and further quantitative analysis.
            pointvalue.Correspondence[ClusterIdx] =  (VoronoiMap->GetPixel(ind));
            MyLabelsOnTrajectory.push_back(VoronoiMap->GetPixel(ind));
            mesh->SetPointData( *pit, pointvalue );
          }
          else
          {
            std::cout << "Point " << point <<  " is outside the space" <<std::endl;
          }

          ++pit;
        }
        // Now compute the miss-matches between the labels on the centers and trajectories:
        int LabelExisted = 0;
        for (unsigned int m=0; m < MyLabels.size(); ++m)
        {
          for (unsigned int n=0; n < MyLabelsOnTrajectory.size(); ++n)
            if (MyLabels.at(m) == MyLabelsOnTrajectory.at(n))
            {
              LabelExisted ++;
              break;
            }
        }


        VariableType AveDist = sumdist/(atrajectory->GetNumberOfPoints());
        VariableType missMatch = (VariableType)(MyLabels.size() - LabelExisted)/MyLabels.size()*AveDist;
        //2nd output
        DissimilarityMeasure[t][ClusterIdx] = (AveDist + missMatch);
        // std::cout << " " << sumdist << " " << atrajectory->GetNumberOfPoints() << " "<< AveDist<<" "<< MyLabels.size()<< " "<<LabelExisted<<" "<< missMatch<<std::endl;
        LabelExisted = 0;

      }
    }
    else
    {
      std::cout<< "Center "<< ClusterIdx <<" is out of space."<< std::endl;
    }

  }

  return DissimilarityMeasure; //NxK
}


VariableType Gamma(VariableType x, VariableType alpha, VariableType beta)
{
  VariableType gammaPdf;

  gammaPdf = 1/(vnl_gamma(alpha)*pow(beta,alpha))*pow(x,(alpha-1))*exp(-x/beta);

  return gammaPdf;
}

Array2DType ComputeLikelihood(const Array2DType &DissimilarityMatrix, ArrayType alpha, ArrayType beta)
{
  unsigned int NumberOfClusters = DissimilarityMatrix.cols();
  unsigned int NumberOfTrajectories = DissimilarityMatrix.rows();

  Array2DType Likelihood;
  Likelihood.SetSize(NumberOfTrajectories,NumberOfClusters);

  for (unsigned int k=0; k<NumberOfClusters; ++k)
  {
    if (beta[k]>0) //valid distrubution
    {
      for (unsigned long int n=0; n<NumberOfTrajectories; ++n)
      {
        VariableType x = DissimilarityMatrix(n,k);
        Likelihood[n][k] = Gamma(x,alpha[k],beta[k]);
      }
    }
    else
    {
      const VariableType z = 0;
      Likelihood.set_column(k,z);
    }
  }
  return Likelihood; //NxK
}

MeshType::Pointer RefineData(const MeshType* mesh, Array2DType &DissimilarityMatrix, Array2DType &Likelihood, Array2DType &Prior, ArrayType MyMinLikelihoodThr, bool havePrior)
{

  MeshType::Pointer RefinedMesh = MeshType::New();
  unsigned long pitNew = 0;
  unsigned long citNew = 0;
  ArrayType arow, peaks;
  Array2DType newDissimilarityMatrix, newLikelihood, newPrior;

  newDissimilarityMatrix.SetSize(DissimilarityMatrix.rows(),DissimilarityMatrix.cols());
  newLikelihood.SetSize(Likelihood.rows(),Likelihood.cols());
  newPrior.SetSize(Prior.rows(),Prior.cols());
  //Store the peak value of likelihood function for each cluster
  peaks.SetSize(Likelihood.cols());
  for (unsigned long int cc=0; cc<Likelihood.cols(); ++cc)
  {
	  peaks[cc] = Likelihood.get_column(cc).max_value();
	  //std::cout << peaks[cc] << std::endl;
  }

  for (unsigned long int n=0; n<Likelihood.rows(); ++n) //go over the trajectories
  {
    arow = Likelihood.get_row(n);
    //To Do: set the threshold differently when havePrior = 1
    bool copycell = 0;
    for (unsigned int m=0; m < arow.Size(); ++m)  //go over the clusters
    if (arow(m)> MyMinLikelihoodThr(m))
      { // To Do: || arow(m)<peaks(m)"
        copycell = 1;
        break;
      }

      if (copycell)
      {
        newDissimilarityMatrix.set_row(citNew, DissimilarityMatrix.get_row(n));
        newLikelihood.set_row(citNew, Likelihood.get_row(n));
        newPrior.set_row(citNew, Prior.get_row(n));


        CellAutoPointer atrajectory, copiedtrajectory;
        copiedtrajectory.TakeOwnership( new PolylineType );

        mesh->GetCell(n, atrajectory);
        PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
        MeshType::PointType point;
        MeshType::PixelType pointvalue;
        CellDataType cellvalue;

        for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
        {
          mesh->GetPoint(*pit, &point);
          mesh->GetPointData( *pit, &pointvalue );
          RefinedMesh->SetPoint(pitNew, point);
          RefinedMesh->SetPointData( pitNew, pointvalue );
          copiedtrajectory->SetPointId(j,pitNew);
          ++pit; ++pitNew;
        }

        mesh->GetCellData(n, &cellvalue);
        RefinedMesh->SetCell(citNew, copiedtrajectory);
        RefinedMesh->SetCellData(citNew, cellvalue);
        citNew++;
      }
  }

  newDissimilarityMatrix = newDissimilarityMatrix.extract(citNew, newDissimilarityMatrix.cols());
  newLikelihood = newLikelihood.extract( citNew, newLikelihood.cols());
  newPrior = newPrior.extract( citNew, newPrior.cols());

  DissimilarityMatrix.SetSize( newDissimilarityMatrix.rows(), newDissimilarityMatrix.cols());
  Likelihood.SetSize( newLikelihood.rows(), newLikelihood.cols());
  Prior.SetSize( newPrior.rows(), newPrior.cols());

  DissimilarityMatrix.copy_in(newDissimilarityMatrix.data_block());
  Likelihood.copy_in(newLikelihood.data_block());
  Prior.copy_in(newPrior.data_block());

  return RefinedMesh;
}
Array2DType ComputePosterior(const Array2DType &Likelihood, const Array2DType &Prior)
{
  Array2DType Posterior;
  Posterior = element_product(Likelihood, Prior);
  //now normalize it:
  if (Likelihood.cols() == 1)
  {
    Posterior /= Posterior.max_value();
  }

  else if (Likelihood.cols() > 1)
  {
    for (unsigned int r = 0; r<Posterior.rows(); ++r)
    {
      VariableType n = Posterior.get_row(r).sum();
      if (n>0)
      {
        Posterior.scale_row(r, 1/n);
      }
    }
  }

  //std::cout<< Posterior << std::endl;
  return Posterior; //NxK
}

void UpdateModelParameters(const Array2DType &DissimilarityMatrix, const Array2DType &Posterior, ArrayType& alpha, ArrayType& beta, Array2DType& W, bool havePrior)
{

  Array2DType A,B;
  A = element_product(Posterior, DissimilarityMatrix);
  B = element_product(Posterior, DissimilarityMatrix.apply(logf));
  unsigned int NumberOfClusters = DissimilarityMatrix.cols();
  ArrayType X, N;
  X.SetSize(NumberOfClusters);
  N.SetSize(NumberOfClusters);

  for (unsigned int k=0; k<NumberOfClusters; ++k)
  {
    N[k] = Posterior.get_column(k).sum();
    if (N[k]>0)
    {
      X[k] = log(A.get_column(k).sum()/N(k))- B.get_column(k).sum()/N(k);
      alpha[k] = (3- X(k) + sqrt(pow(X(k)-3,2)+ 24*X(k)))/(12*X(k));
      beta[k] = A.get_column(k).sum()/(alpha(k)*N(k));
      if (!havePrior)
      { //update the mixing weights
        W.set_column(k, N[k]/Posterior.rows());
      }
    }
    else
    { //null cluster - no distribution
      alpha[k] = 0;
      beta[k] = 0;
      const VariableType z = 0;
      W.set_column(k,z);
    }
  }
  // std::cout << W <<std::endl;
}
void AssignClusterLabels(MeshType* mesh, const Array2DType &Posterior)
{
  if (Posterior.rows()!= mesh->GetNumberOfCells())
  {
    std::cerr<< "There is a miss-match between the number of trajectories and the membership probabilities to be assigned" <<std::endl;
  }
  CellDataType cellvalue;
  ArrayType arow;
  for (unsigned int t=0; t<Posterior.rows(); ++t)
  {
    VariableType my_max = -1;
    int my_max_idx;
    arow = Posterior.get_row(t);
    // find the maximum membership probability for t'th trajectory
    for (unsigned int m=0; m<arow.Size(); ++m)
    {
      if (arow(m)>my_max)
      {
        my_max = arow(m);
        my_max_idx = m;
      }
    }
    CellAutoPointer atrajectory;
    mesh->GetCell(t, atrajectory);
    mesh->GetCellData(t, &cellvalue );
    cellvalue.ClusterLabel = my_max_idx+1;    //start the Cluster Labels (Ids) from 1
    cellvalue.membershipProbability.SetSize(arow.Size());
    cellvalue.membershipProbability = arow;
    mesh->SetCellData(t, cellvalue );
  }
}

MeshType::Pointer UpdateCenters(MeshType* mesh, MeshType* mesh_centers, const Array2DType &Posterior, ImageType* refImage, VariableType MinPosterior, VariableType MaxDist)
{
  MeshType::Pointer mesh_newcenters=MeshType::New();
  ArrayType post;
  unsigned int NumberOfClusters=mesh_centers->GetNumberOfCells();
  long int cit =0;

  for (unsigned int k=0; k<NumberOfClusters; ++k)
  {
    post = Posterior.get_column(k);
    std::vector<long int> sigIDs;
    for (unsigned long int p=0; p<post.Size(); ++p)
    {
      if (post(p)>MinPosterior)
      {
        sigIDs.push_back(p);
      }
    }

    CellAutoPointer Centerline, new_center;
    new_center.TakeOwnership( new PolylineType );

    mesh_centers->GetCell(k,Centerline);
    PolylineType::PointIdIterator mcit = Centerline->PointIdsBegin();
    ImageType::IndexType ind;
    MeshType::PointType point, last_mean_point;
    VariableType distBetweenSuccessivePoints =0.5;


    int MyLabel;
    unsigned int s=0;
    for (unsigned int c=0; c<Centerline->GetNumberOfPoints(); ++c)
    {
      mesh_centers->GetPoint(*mcit, &point);
      refImage->TransformPhysicalPointToIndex(point, ind);
      MyLabel = refImage->ComputeOffset(ind);

      MeshType::PixelType tpointvalue;
      std::vector<MeshType::PointType> pntStack;
      std::vector<VariableType> postStack;

      MeshType::PointType tpoint, mean_point,sum_points, closest_point;
      VariableType dist,closest_point_post;
      pntStack.clear();
      postStack.clear();

      if (sigIDs.size()>0)
      {
        //update the center by taking the average of trajectories

        for (unsigned long int t=0; t<sigIDs.size(); ++t)
        {
          CellAutoPointer atrajectory;
          mesh->GetCell(sigIDs.at(t),atrajectory);
          PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();

          VariableType MinDist = itk::NumericTraits<VariableType>::max();
          for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
          {
            mesh->GetPoint(*pit, &tpoint);
            mesh->GetPointData( *pit, &tpointvalue );
            //find the points on a single trajectory that correponde to the current point on the center
            if (tpointvalue.Correspondence(k)==MyLabel)
            {
              dist = tpoint.EuclideanDistanceTo(point);
              if (dist<MinDist)
              {
                MinDist = dist;
                closest_point = tpoint;
                closest_point_post = post(sigIDs.at(t));
              }
            }
            pit++;
          }
          if (MinDist<MaxDist)
          {
            pntStack.push_back(closest_point);
            postStack.push_back(closest_point_post);
          }

        }
        // if (pntStack.size()<3)
        //   std::cout<<"Point "<< c <<" on the new center is obtained by aveaging less than 3 points." << std::endl;
        sum_points.Fill(0);
        VariableType SumPost = 0;
        for ( unsigned int m=0; m<pntStack.size(); ++m)
        {
          MeshType::PointType temp = pntStack[m];
          sum_points[0] += temp[0]*postStack[m];
          sum_points[1] += temp[1]*postStack[m];
          sum_points[2] += temp[2]*postStack[m];
          SumPost +=postStack[m];
        }
        if (SumPost>0)
        {
          mean_point[0] = sum_points[0]/SumPost;
          mean_point[1] = sum_points[1]/SumPost;
          mean_point[2] = sum_points[2]/SumPost;

          //compute the distance between the current mean point and the previous one:
          if (c>0 && cit>0) //not at the begining of the centerline
          {
            distBetweenSuccessivePoints = mean_point.EuclideanDistanceTo(last_mean_point);
          }
          if (distBetweenSuccessivePoints>=0.5)     /////// parameter: dist between the samples!
          {
            mesh_newcenters->SetPoint(cit,mean_point);
            new_center->SetPointId(s,cit);
            last_mean_point = mean_point;
            ++cit; ++s;
          }
        }
        else
        {
          //std::cout<<"A point on the center is not being updated!" << std::endl;
        }
        ++mcit;
      }

      else
      {//juts copy the points from the center to newcenters
        mesh_newcenters->SetPoint(cit,point);
        new_center->SetPointId(s,cit);
        ++cit; ++s; ++mcit;
      }
    }

    mesh_newcenters->SetCell(k,new_center);
  }

  return mesh_newcenters;
}

CurveType SmoothCurve(CurveType Curve)
{
  CurveType SmoothedCurve;
  int NumberOfPoints = Curve.rows();
  SmoothedCurve.set_size(NumberOfPoints,3);
  int window = 5;  //radius
  for (int j=0; j<NumberOfPoints; ++j)
  {
    CurvePointType sumPoints;
    sumPoints.set_size(3);
    sumPoints.fill(0);
    int el =0,low,up;
    if  ((j-window)<0)
    {
      low=0;
    }
    else
    {
      low =(j-window);
    }

    if  (NumberOfPoints<(j+window)) up=NumberOfPoints; else up =(j+window);
    for (int i=low; i< up; ++i)
    {
      sumPoints += Curve.get_row(i); el++;
    }
    SmoothedCurve.set_row(j, sumPoints/el);
  }
  //std::cout << SmoothedCurve <<std::endl;
  return SmoothedCurve;
}

MeshType::Pointer SmoothMesh(MeshType* mesh)
{
  MeshType::Pointer smoothedMesh = MeshType::New();
  unsigned int NumberOfCells = mesh->GetNumberOfCells();
  long unsigned int newpid = 0;
  for (unsigned int k=0; k<NumberOfCells; ++k)
  {
    CellAutoPointer aCell,newCell;
    newCell.TakeOwnership( new PolylineType );

    mesh->GetCell(k,aCell);

    MeshType::PointType point;
    CurveType MyCurve, SmoothedCurve;
    MyCurve.SetSize(aCell->GetNumberOfPoints(),3);

    PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
    for (unsigned int j=0; j < aCell->GetNumberOfPoints(); ++j)
    {
      mesh->GetPoint(*pit, &point);
      MyCurve.set_row(j, point.GetVnlVector());
      pit++;
    }

    SmoothedCurve = SmoothCurve(MyCurve);
    //Put the new curve in the mesh:
    for (unsigned int j=0; j < SmoothedCurve.rows(); ++j)
    {
      MeshType::PointType mpoint;
      mpoint[0] = SmoothedCurve(j,0);
      mpoint[1] = SmoothedCurve(j,1);
      mpoint[2] = SmoothedCurve(j,2);
      smoothedMesh->SetPoint(newpid, mpoint);
      newCell->SetPointId(j,newpid);
      newpid++;
    }
    smoothedMesh->SetCell(k,newCell);
  }
  return smoothedMesh;
}

VariableType diffCurve(CurveType MyCurve1, CurveType MyCurve2) //simplest implementation
{
  VariableType dist = 0;
  CurvePointType p1,p2;
  for (unsigned int l = 0; l< MyCurve2.rows(); ++l)
  {
    p2 = MyCurve2.get_row(l);
    p1 = MyCurve1.get_row(l);
    dist+=(p1 - p2).two_norm();
  }
  return dist/MyCurve2.rows();
}

ArrayType diffMeshes(const MeshType* mesh1, const MeshType* mesh2)
{
  ArrayType dist;

  unsigned int NumberOfCells1 = mesh1->GetNumberOfCells();
  unsigned int NumberOfCells2 = mesh2->GetNumberOfCells();
  dist.SetSize(NumberOfCells1); dist.fill(10000);
  if (NumberOfCells1 != NumberOfCells2)
  {
    std::cout << "Number of cells don't match between the given meshes" <<std::endl;
  }
  else
  {
    dist.fill(0);
    //long unsigned int newpid = 0;
    for (unsigned int k=0; k<NumberOfCells1; ++k)
    {
      CellAutoPointer aCell1,aCell2;
      //newCell.TakeOwnership( new PolylineType );

      mesh1->GetCell(k,aCell1);
      mesh2->GetCell(k,aCell2);

      MeshType::PointType point;
      CurveType MyCurve1, MyCurve2;
      MyCurve1.SetSize(aCell1->GetNumberOfPoints(),3);
      MyCurve2.SetSize(aCell2->GetNumberOfPoints(),3);

      PolylineType::PointIdIterator pit1 = aCell1->PointIdsBegin();
      for (unsigned int j=0; j < aCell1->GetNumberOfPoints(); ++j)
      {
        mesh1->GetPoint(*pit1, &point);
        MyCurve1.set_row(j, point.GetVnlVector());
        pit1++;
      }

      PolylineType::PointIdIterator pit2 = aCell2->PointIdsBegin();
      for (unsigned int i=0; i < aCell2->GetNumberOfPoints(); ++i)
      {
        mesh2->GetPoint(*pit2, &point);
        MyCurve2.set_row(i, point.GetVnlVector());
        pit2++;
      }


      dist[k] = diffCurve(MyCurve1, MyCurve2);

    }
  }
  return dist;

}



void SetInitialValue(const Array2DType &DissimilarityMatrix, ArrayType &beta)
{
  VariableType MaxErr = 10;
  for (unsigned int k = 0; k<DissimilarityMatrix.cols(); ++k)
  {
    VariableType sumErr = 0; long int n =0;
    for (unsigned long int t = 0; t<DissimilarityMatrix.rows(); ++t)
      if (DissimilarityMatrix(t,k)<MaxErr)
      {
        sumErr+=DissimilarityMatrix(t,k);
        n++;
      }
      beta[k] = sumErr/n;
      //std::cout<< beta[k] << std::endl;
      //beta[k] = 5;
  }
}

ArrayType AdjustThreshold(VariableType factor, ArrayType alpha, ArrayType beta)
{
  ArrayType MyMinLikelihoodThr;
  MyMinLikelihoodThr.set_size(alpha.size());
  for (unsigned int k=0; k<alpha.size(); ++k)
  {
    if (beta[k]>0) //valid distribution
    {
      VariableType GammaMode = (alpha[k]>=1)?((alpha[k]-1)*beta[k]):0.1;
      MyMinLikelihoodThr[k] = factor*Gamma(GammaMode,alpha[k],beta[k]);
    }
    else
    {
      MyMinLikelihoodThr[k] = itk::NumericTraits<VariableType>::max();
    }
  }
  return MyMinLikelihoodThr;
}

void ComputeScalarMeasures(MeshType* Trajectories)
{
  for (unsigned long int t=0; t<Trajectories->GetNumberOfPoints(); ++t)
  {
    MeshType::PixelType pointvalue;
    itk::FixedArray<double, 9 > tensor9;
    TensorPixelType tensor6;

    Trajectories->GetPointData(t, &pointvalue );
    tensor9 = pointvalue.Tensor;
    tensor6[0]=tensor9[0];
    tensor6[1]=tensor9[1];
    tensor6[2]=tensor9[2];
    tensor6[3]=tensor9[4];
    tensor6[4]=tensor9[5];
    tensor6[5]=tensor9[8];
    EigenValuesArrayType eigenvals;
    tensor6.ComputeEigenValues(eigenvals);
    pointvalue.EigenValues = eigenvals;
    pointvalue.FA = tensor6.GetFractionalAnisotropy();

    if (pointvalue.FA>=1)
    {
      pointvalue.FA = 0.5; // for now
      std::cerr << " Reached at a point with FA greater than 1" <<std::endl;
    }
    Trajectories->SetPointData(t, pointvalue );

  }
}

Array3DType BuildFeatureMatrix(const MeshType* cluster, const MeshType* center, int clusterId, ImageType* refImage)
{
  Array2DType fMatrix1,fMatrix2,fMatrix3,fMatrix4;  // NxS (Number of Trajectories x Number of Samples on the Center)
  VariableType nanVal = -1;
  unsigned long int numberOfTrajectories = cluster->GetNumberOfCells();
  unsigned int numberOfSamples = center->GetNumberOfPoints();
  fMatrix1.set_size(numberOfTrajectories,numberOfSamples);
  fMatrix2.set_size(numberOfTrajectories,numberOfSamples);
  fMatrix3.set_size(numberOfTrajectories,numberOfSamples);
  fMatrix4.set_size(numberOfTrajectories,numberOfSamples);

  ImageType::IndexType ind;
  MeshType::PointType point;
  int MyLabel, PrevLabel=0;
  for (unsigned int s=0; s<numberOfSamples; ++s)
  {
    center->GetPoint(s, &point);
    if (refImage->TransformPhysicalPointToIndex(point, ind))
    {
      MyLabel = refImage->ComputeOffset(ind);
    }
    else
    {
      std::cout<< "point is out of space"<< std::endl;
    }
    if (MyLabel!=PrevLabel) // Because of quantization, labels of successive samples could be the same.
    {
      PrevLabel = MyLabel;
      //go over trajectories
      for (unsigned long int t=0; t<numberOfTrajectories; ++t)
      {
        CellAutoPointer atrajectory;
        cluster->GetCell(t,atrajectory);
        PolylineType::PointIdIterator pit = atrajectory->PointIdsBegin();
        MeshType::PointType tpoint;
        MeshType::PixelType tpointvalue;

        double sumFeature1 = 0;
        double sumFeature2 = 0;
        double sumFeature3 = 0;
        double sumFeature4 = 0;
        int n=0;
        for (unsigned int j=0; j < atrajectory->GetNumberOfPoints(); ++j)
        {
          cluster->GetPoint(*pit, &tpoint);
          cluster->GetPointData( *pit, &tpointvalue );
          if (tpointvalue.Correspondence(clusterId)==MyLabel)
          {
            sumFeature1 += tpointvalue.FA;
            sumFeature2 += (tpointvalue.EigenValues[0]+tpointvalue.EigenValues[1]+tpointvalue.EigenValues[2])/3;
            sumFeature3 += (tpointvalue.EigenValues[0]+tpointvalue.EigenValues[1])/2;
            sumFeature4 += tpointvalue.EigenValues[2];
            n++;
          }
          pit++;
        }
        if (sumFeature1>0)
        {
          fMatrix1[t][s] = sumFeature1/n;
          fMatrix2[t][s] = sumFeature2/n;
          fMatrix3[t][s] = sumFeature3/n;
          fMatrix4[t][s] = sumFeature4/n;

        }
        else
        {
          fMatrix1[t][s] = nanVal;
          fMatrix2[t][s] = nanVal;
          fMatrix3[t][s] = nanVal;
          fMatrix4[t][s] = nanVal;
        }
      }
    }
    else
    {
      for (unsigned long int t=0; t<numberOfTrajectories; ++t)
      {
        fMatrix1[t][s] = fMatrix1(t,s-1);
        fMatrix2[t][s] = fMatrix2(t,s-1);
        fMatrix3[t][s] = fMatrix3(t,s-1);
        fMatrix4[t][s] = fMatrix4(t,s-1);
      }
    }
  }
  Array3DType allFeatures;
  allFeatures.push_back(fMatrix1);
  allFeatures.push_back(fMatrix2);
  allFeatures.push_back(fMatrix3);
  allFeatures.push_back(fMatrix4);

  return allFeatures;
}

MeshType::Pointer getTrajectories(MeshType* mesh, std::vector<unsigned long int> CellIDs)
{
  MeshType::Pointer selectedCells=MeshType::New();
  CellAutoPointer aCell, MyCell;
  long int myid =0;

  for (unsigned int c=0; c<CellIDs.size(); ++c)
  {
    MyCell.TakeOwnership( new PolylineType );
    MeshType::CellIdentifier CellID = CellIDs.at(c);
    mesh->GetCell(CellID, aCell);
    CellDataType cellvalue;
    mesh->GetCellData(CellID, &cellvalue);
    PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
    MeshType::PointType point;
    MeshType::PixelType pointvalue;
    for (unsigned int j=0; j < aCell->GetNumberOfPoints(); j++)
    {
      mesh->GetPoint(*pit, &point);
      mesh->GetPointData(*pit, &pointvalue);
      selectedCells->SetPoint(myid, point );
      selectedCells->SetPointData(myid, pointvalue );
      MyCell->SetPointId(j,myid);
      pit++; myid++;
    }
    selectedCells->SetCell(c, MyCell);
    selectedCells->SetCellData(c, cellvalue);

  }

  return selectedCells;
}

MeshType::Pointer getCluster(MeshType* Trajectories,int k)
{
  MeshType::Pointer cluster;
  std::vector<unsigned long int> CellIds;
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellAutoPointer atrajectory;
    CellDataType cellvalue;
    Trajectories->GetCell(t, atrajectory);
    Trajectories->GetCellData(t, &cellvalue );
    if (cellvalue.ClusterLabel == k+1)
    {
      CellIds.push_back(t);
    }
  }
  cluster = getTrajectories(Trajectories, CellIds);
  return cluster;
}

MeshType::Pointer getCluster(MeshType* Trajectories, std::string caseName)
{
  MeshType::Pointer cluster;
  std::vector<unsigned long int> CellIds;
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellAutoPointer atrajectory;
    CellDataType cellvalue;
    Trajectories->GetCell(t, atrajectory);
    Trajectories->GetCellData(t, &cellvalue );
    if (cellvalue.CaseName.compare(caseName)==0 )
    {
      CellIds.push_back(t);
    }
  }
  cluster = getTrajectories(Trajectories, CellIds);
  return cluster;
}

Array2DType getClusterPosterior(Array2DType Posterior, MeshType* Trajectories,int k)
{
  Array2DType post;
  std::vector<long int> CellIds;
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellAutoPointer atrajectory;
    CellDataType cellvalue;
    Trajectories->GetCell(t, atrajectory);
    Trajectories->GetCellData(t, &cellvalue );
    if (cellvalue.ClusterLabel == k+1)
    {
      CellIds.push_back(t);
    }
  }
  post.set_size(CellIds.size(),1);
  for (unsigned long int c=0; c<CellIds.size(); ++c)
  {
    post.set_row(c,Posterior.get(CellIds[c],k)); //Posterior.get_row(CellIds[c])
  }
  return post;
}

std::vector<std::string> getClusterSubjectNames(MeshType* Trajectories)
{
	std::vector<std::string> subjectNames;
	for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellDataType cellvalue;
    Trajectories->GetCellData(t, &cellvalue );
    subjectNames.push_back(cellvalue.CaseName);
  }

  return subjectNames;
}


ImageType::Pointer getSubSpace(const MeshType* Trajectories, std::vector<double> spacing)
{
  ImageType::Pointer subSpace = ImageType::New();
  MeshType::BoundingBoxType const *dataBoundingBox = Trajectories->GetBoundingBox();
  dataBoundingBox->ComputeBoundingBox();
  MeshType::BoundingBoxType::BoundsArrayType dataBounds = dataBoundingBox->GetBounds();

  MeshType::PointType p1,p2,p;

  p1[0] = dataBounds[0];
  p1[1] = dataBounds[2];
  p1[2] = dataBounds[4];
  p2[0] = dataBounds[1];
  p2[1] = dataBounds[3];
  p2[2] = dataBounds[5];

  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  ImageType::SizeType size;
  size[0] = ceil((p2[0]-p1[0]+5)/spacing[0]);
  size[1] = ceil((p2[1]-p1[1]+5)/spacing[1]);
  size[2] = ceil((p2[2]-p1[2]+5)/spacing[2]);

  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );

  subSpace->SetRegions(desiredRegion);
  subSpace->Allocate();

  ImageType::SpacingType imSpacing;
  imSpacing[0] = spacing[0];
  imSpacing[1] = spacing[1];
  imSpacing[2] = spacing[2];


  subSpace->SetSpacing(imSpacing);
  subSpace->SetOrigin(p1);

  return subSpace;
}

void  fillPriorInfo(Array2DType &Prior, MeshType* Trajectories)
{
  for (unsigned long int t=0; t<Trajectories->GetNumberOfCells(); ++t)
  {
    CellDataType cellvalue;
    Trajectories->GetCellData(t, &cellvalue);
    Prior.set_row(t, cellvalue.atlasPriors);
  }
}

ArrayType meanMat(const Array2DType &X, int nanVal=0)
//take the column-wise mean of the matrix X, ignoring the zero elements.
{
  ArrayType mX;
  mX.SetSize(X.cols());
  ArrayType aCol;
  for (unsigned int c = 0; c<X.cols(); c++)
  {

    aCol = X.get_column(c);

    // if there is no NAN in the colume:   (here we assume that the nanVal is 0 or a negative number.
    if (aCol.min_value() != nanVal)
    {
      mX(c) = aCol.mean();
    }
    else
    {
      double s = 0;
      unsigned int n = 0;
      for (unsigned int i =0; i<aCol.Size(); i++)
      {
        if (aCol(i)!= nanVal)
        {
          s+=aCol(i);
          n++;
        }
      }
      if (n!=0)
      {
        mX(c) = s/n;
      }
      else
      {
        //
        mX(c)=0;
        std::cout << "NaN column at " << c << "!" <<std::endl;
      }
    }

  }

  return mX;
}

ArrayType stdMat(const Array2DType &X, int nanVal=0)
//take the column-wise std of the matrix X, ignoring the nonVal elements.
{
  ArrayType mX;
  mX.SetSize(X.cols());
  ArrayType aCol;
  for (unsigned int c = 0; c<X.cols(); c++)
  {

    aCol = X.get_column(c);

    // if there is no NAN in the colume:   (here we assume that the nanVal is 0 or a negative number.
    if (aCol.min_value() != nanVal)
    {
		mX(c) = (aCol - aCol.mean()).rms();
    }
    else
    {
      double s = 0;
      unsigned int n = 0;
      for (unsigned int i =0; i<aCol.Size(); i++)
      {
        if (aCol(i)!= nanVal)
        {
          s+=aCol(i);
          n++;
        }
      }
      if (n!=0)
      {
         double meanVal =  s/n;
		 double ms = 0;
		 for (unsigned int i =0; i<aCol.Size(); i++)
			{
			if (aCol(i)!= nanVal)
			{
				ms+=(aCol(i)-meanVal)*(aCol(i)-meanVal);

			}
			}


		mX(c) = sqrt(ms/n);

      }
      else
      {
        //
        mX(c)=0;
        std::cout << "NaN column at " << c << "!" <<std::endl;
      }
    }

  }

  return mX;
}


ArrayType meanMat(Array2DType X, Array2DType P, int nanVal=0)
//take the column-wise 'weighted mean' of the matrix X, ignoring the zero elements.
{
  X = element_product(X, P);   //x = x.p(x)
  ArrayType mX;
  mX.SetSize(X.cols());
  ArrayType aCol;
  for (unsigned int c = 0; c<X.cols(); c++)
  {

    aCol = X.get_column(c);
    VariableType s;

    if (aCol.min_value() != nanVal)
    {
      s = aCol.sum();
      mX(c) = s/P.get_column(c).sum();
    }
    else
    {
      VariableType n = 0;
      for (unsigned int i =0; i<aCol.Size(); i++)
      {
        if (aCol(i) != nanVal)
        {
          n+= P.get_column(c).get(i);
          s+= aCol(i);
        }

      }
      if (n!=0)
      {
        mX(c) = s/n;
      }
      else
      {
        mX(c)=0;
        std::cout << "NaN column at " << c << "!" <<std::endl;
      }
    }

  }

  return mX;
}

void  AddPointScalarToACell(MeshType* mesh, MeshType::CellIdentifier CellID, ArrayType f)
{
  CellAutoPointer aCell;
  mesh->GetCell(CellID, aCell);
  PolylineType::PointIdIterator pit = aCell->PointIdsBegin();
  MeshType::PixelType pointvalue;
  for (unsigned int j=0; j < aCell->GetNumberOfPoints(); j++)
  {
    mesh->GetPointData(*pit, &pointvalue);
    pointvalue.FA = f[j];
    mesh->SetPointData(*pit, pointvalue );
    pit++;
  }

}

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


int main(int argc, char* argv[])
{

  PARSE_ARGS;

  MeshType::Pointer    Trajectories, Centers, atlasCenters;
  MeshType::Pointer    oldCenters = MeshType::New();

  populationStudy=0;
  useAtlas = 0;
  PerformQuantitativeAnalysis = 1;
  CopyFieldType copyField = {0,0,1,1,0};

  // Get the input trajectories

  std::vector<std::string> allfilenames;
  if (populationStudy)
  {
	copyField.CaseName = 1;
	itksys::Glob allfiles;
	std::string globPath = TractsDir+"/*.vt*"; //vtk or vtp
	allfiles.FindFiles(globPath);
	allfilenames = allfiles.GetFiles();
	Trajectories = ReadVTKfiles(allfilenames);
  }
  else if (!trajectoriesFilename.empty())
  {
	  copyField.CaseName = 0;
      Trajectories = ReadVTKfile(trajectoriesFilename.c_str());
  }
  else
  {
	std::cerr << "No input was given as a collection of trajectories to be clustered" << std::endl;
	return -1;
  }

  // Get the intialcenter(s)
  if (useAtlas)
  {
    std::string path;
    if (AtlasDir[0] == '.')
    {
      path = itksys::SystemTools::GetFilenamePath(argv[0]);
      path = path + "/" + AtlasDir;
    }
    char atlasfileName1[250];
    sprintf(atlasfileName1, "%s/atlasCenters.vtp",path.c_str());
	atlasCenters = ReadVTKfile(atlasfileName1);

    //Read FA map of the atlas
    char atlasfileName2[250];
	sprintf(atlasfileName2, "%s/atlasFAMap.nhdr",path.c_str());
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(atlasfileName2);
    reader->Update();
    ImageType::Pointer atlasFAVolume = reader->GetOutput();
    //Read FA map of the case
    ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName(subjectFAFilename.c_str());
    reader2->Update();
    ImageType::Pointer caseFAVolume = reader2->GetOutput();
    //Read which bundles are selected to be clustered:
    std::vector<unsigned long int> atlasCellIDs;
    if (splenium)
    {atlasCellIDs.push_back(0);}
    if (genu)
    {atlasCellIDs.push_back(2);}
    if (cingulumR)
    {atlasCellIDs.push_back(3);}
    if (cingulumL)
    {atlasCellIDs.push_back(4);}

    if (atlasCellIDs.size()>0)
    {
      //Do affine registration
      TransformType::Pointer transform;
      transform = doAffineRegistration(caseFAVolume, atlasFAVolume);
      //Select and transfer the centers
       Centers = applyTransform(atlasCenters, transform, atlasCellIDs );
       CopyFieldType copyField = {0,0,0,0};
       WriteVTKfile(Centers,transformedCentersFilename.c_str(),copyField);
    }
    else
    {
      std::cerr << "User should select at least one of the bundles in the atlas' list to be clustered! " << std::endl;
      return -1;
    }
  }
  else if (!centersFilename.empty())
  {
    Centers = ReadVTKfile(centersFilename.c_str());
  }

  else
    {
    std::cerr << "No initial center was given" << std::endl;
    return -1;
    }

  //Write initial centers:

  copyField.ClusterLabel = 0;

  //WriteVTKfile(Centers,initCentersFilename.c_str(),copyField);

  copyField.ClusterLabel = 1;
  VariableType MinPost = (VariableType) 1/(Centers->GetNumberOfCells());
  VariableType MinLike = 0.1*MinLikelihoodThr;   // 5->0.5 ; 1 ->0.1

  ArrayType alpha, beta, MyMinLikelihoodThr;
  Array2DType DissimilarityMatrix, Likelihood, Prior, Posterior; //NxK

  //EM Initialization
  alpha.SetSize(Centers->GetNumberOfCells()); alpha.fill(1);
  beta.SetSize(Centers->GetNumberOfCells()); beta.fill(5);
  Prior.SetSize(Trajectories->GetNumberOfCells(),Centers->GetNumberOfCells());
  bool havePrior = 0;
  if (havePrior)
  {
    fillPriorInfo(Prior, Trajectories);
  }
  else //in absence of an atlas
  {
    VariableType initp = 1.0 / (Centers->GetNumberOfCells());
    Prior.fill(initp);
  }

  // set the space to the limits of input trajectories
  ImageType::Pointer subSpace;
  subSpace = getSubSpace(Trajectories, spacing);

  ///// START /////

  bool debug = 0;
  ArrayType dd;
  for (int i=0; i<maxNumberOfIterations; ++i)
  {
    std::cout<< "Iteration  " << i+1 << std::endl;
    DissimilarityMatrix = ComputeDissimilarity(Trajectories, Centers, subSpace);
    if (debug)
    {
      std::cout<< DissimilarityMatrix << std::endl;
      std::cout<< "beta = " << beta << std::endl;
    }

    //EM Block Begins
    Likelihood = ComputeLikelihood(DissimilarityMatrix, alpha, beta);
    MeshType::Pointer RefinedTrajectories;
    MyMinLikelihoodThr = AdjustThreshold(MinLike, alpha, beta);

    if (debug)
    {
      std::cout<< Likelihood << std::endl;
      std::cout<< "MyMinLikelihoodThr = " << MyMinLikelihoodThr << std::endl;
    }

    RefinedTrajectories = RefineData(Trajectories,DissimilarityMatrix,Likelihood,Prior,MyMinLikelihoodThr,havePrior);
    if (RefinedTrajectories->GetNumberOfCells()<1)
    {
      std::cerr<< "The current setting of data/parameters have resulted in zero clustered trajectories"<<std::endl;
      return 0; //EXIT_FAILURE;
    }
    Posterior = ComputePosterior(Likelihood,Prior);
    UpdateModelParameters(DissimilarityMatrix,Posterior,alpha,beta,Prior,havePrior);
    //EM Block Ends


    if (debug)
    {
      WriteVTKfile(RefinedTrajectories, "RefinedTraj.vtp", copyField);
      std::cout<< Posterior << std::endl;
      std::cout<< "alpha = " << alpha << std::endl;
      std::cout<< "beta = " << beta << std::endl;

      char fileName[250];
      sprintf(fileName, "%s/posterior%d.csv", OutputDirectory.c_str(),i+1);
      WriteCSVfile(fileName, Posterior);
    }

    //Update centers:
    MeshType::Pointer NewCenters;
    NewCenters = UpdateCenters(RefinedTrajectories, Centers, Posterior, subSpace, MinPost, MaxDist);

    //Smooth centers:
    MeshType::Pointer SmoothedCenters;
    SmoothedCenters = SmoothMesh(NewCenters);

    if (debug)
    {
      char tempcenterName[250];
      sprintf(tempcenterName, "%s/centers_iter%d.vtp", OutputDirectory.c_str(),i+1);
      copyField.FA = 0; copyField.Tensor = 0;
      WriteVTKfile(NewCenters,tempcenterName,copyField);
    }



    oldCenters = Centers;
    Centers = SmoothedCenters;
    Trajectories = RefinedTrajectories;
    if (i>1)
    {
    dd = diffMeshes(oldCenters, Centers);
	  //std::cout<< "Difference between new centers and old ones is "<< dd.max_value() <<std::endl;
    if (dd.max_value()<5) break;
    }
  }

  AssignClusterLabels(Trajectories,Posterior);
  copyField.FA = 0; copyField.Tensor = 1;
  WriteVTKfile(Trajectories, outputClustersFilename.c_str(),copyField);
  //Done with clustering.


  //////////////////////////////////////////////////////////////////////
  //Start Quantitative Analysis:
  //////////////////////////////////////////////////////////////////////
  if (PerformQuantitativeAnalysis)
  {
    //Compute and add diffusion scalar measures to each point on the mesh:
    ComputeScalarMeasures(Trajectories);

    //Generate separate mesh for each cluster -> cluster + center:
    std::vector <MeshType::Pointer> cluster, center, centerWithData;
    std::vector<unsigned long int> cellId;
    std::vector <Array3DType> clusterFeatures;
    std::stringstream currentClusterId;
    Array3DType posts;

    for (unsigned int k=0; k<Centers->GetNumberOfCells(); ++k)
    {
      currentClusterId<<k+1;

      //Separate cluster k'th
      cluster.push_back(getCluster(Trajectories,k));

	  //Separate center  k'th
      cellId.clear(); cellId.push_back(k);
      center.push_back(getTrajectories(oldCenters,cellId));

	  //Separate posterior probabilities
      posts.push_back(getClusterPosterior(Posterior,Trajectories,k));

	  //Compute the feature matrices
      clusterFeatures.push_back(BuildFeatureMatrix(cluster[k],center[k],k, subSpace));
      if (clusterFeatures[k].at(0).rows()>0) //not empty
	  {
	    //Now compute the mean FA and assign it to the pointvalue.FA of each point on the center
		ArrayType meanFA; //, stdFA;
		meanFA = meanMat(clusterFeatures[k].at(0),-1);                          //TO Do: compute the weighted average
	   //Add point data to the cell with the cell ID of cellId in the oldCenters mesh to get visualized with the
	   // mean FA after loading in Slicer:
		AddPointScalarToACell(oldCenters,k, meanFA );//oldCenters gets updated.
	  }

      centerWithData.push_back(getTrajectories(oldCenters,cellId));

	  //write out the posterior for each cluster
	  std::string filename;
	  if (posts[k].rows()>0)
	  {
		filename = OutputDirectory+"/cluster"+ currentClusterId.str() + "_posterior.csv";
		WriteCSVfile(filename, posts[k]);
		//now write out the feature matrices to files:
		filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_FA.csv";
		WriteCSVfile(filename,clusterFeatures[k].at(0));

		filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_MD.csv";
	    WriteCSVfile(filename,clusterFeatures[k].at(1));

	    filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_PerDiff.csv";
	    WriteCSVfile(filename,clusterFeatures[k].at(2));

	    filename = OutputDirectory + "/cluster"+ currentClusterId.str()+"_ParDiff.csv";
	    WriteCSVfile(filename,clusterFeatures[k].at(3));

       // Write individual files for each cluster and its center.
        filename = OutputDirectory+"/center" + currentClusterId.str()+".vtp";
		copyField.FA = 1; copyField.Tensor = 0;
		WriteVTKfile(centerWithData[k],filename,copyField);

		filename = OutputDirectory+"/cluster" + currentClusterId.str()+".vtp";
		copyField.FA = 0; copyField.Tensor = 1;
		WriteVTKfile(cluster[k], filename,copyField);
	  }

	  std::vector<std::string> subjectNames;
      std::vector<MeshType::Pointer> subClusters;

	  if (populationStudy)
	  {

		  std::vector<std::string> subjectNamesInCluster =  getClusterSubjectNames(cluster[k]);
		  std::string myfilename = OutputDirectory + "/cluster" + currentClusterId.str()+ "_subjectNames.txt";
		  ofstream myfile;
		  std::cout<< "Writing out "<< myfilename.c_str() << "..." << std::endl;
		  myfile.open (myfilename.c_str());
		  for (unsigned int m =0; m<subjectNamesInCluster.size(); m++)
		  {
			  myfile << subjectNamesInCluster.at(m);
			  myfile << std::endl;
		  }
	      myfile.close();

		  for (unsigned int sn=0; sn<allfilenames.size(); sn++)
		  {
			  std::string filename = itksys::SystemTools::GetFilenameName(allfilenames.at(sn));
			  std::string subjectName = filename.substr(0,13);
			  std::string outputfilename = OutputDirectory + "/cluster" + currentClusterId.str() + "_" + subjectName + ".vtp";

			  //Generate separate meshes for each subject in the population
		      subClusters.push_back(getCluster(cluster[k], allfilenames.at(sn)));
              WriteVTKfile(subClusters[sn], outputfilename.c_str() ,copyField);

		  }

	  }

    }
		copyField.FA = 1;
    }
  else
  {
	    copyField.FA = 0;
  }

  //write centers with new point data from the quantitative analysis part.
  copyField.Tensor = 0;
  WriteVTKfile(oldCenters, outputCentersFilename.c_str(),copyField);

  std::cout << "Done." << std::endl;
  return EXIT_SUCCESS;
}
