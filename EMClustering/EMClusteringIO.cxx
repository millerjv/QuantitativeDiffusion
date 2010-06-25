
#include <EMClusteringIO.h>


vtkPolyData* itk2vtkPolydata(MeshType* mesh, CopyFieldType copyField)
{
  vtkPolyData* polydata = vtkPolyData::New();
  //CopyItkMesh2VtkPolyData(mesh, polydata, copyField);
  //void CopyItkMesh2VtkPolyData(MeshType* mesh, vtkPolyData* polydata, CopyFieldType copyField)

  // Convert the itk mesh to vtk polydata:
  //{
    unsigned int numPoints = mesh->GetNumberOfPoints();
    unsigned int numCells  = mesh->GetNumberOfCells();

    vtkPoints* vpoints = vtkPoints::New();
    vpoints->SetNumberOfPoints(numPoints);

    vtkDoubleArray* tensors = vtkDoubleArray::New();
    tensors->SetNumberOfComponents(9);
    tensors->SetNumberOfTuples(numPoints);

    vtkLongArray* correspondences = vtkLongArray::New();
    correspondences->SetNumberOfComponents(1);
    correspondences->SetNumberOfTuples(numPoints);
    correspondences->SetName("Correspondence");

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
    itk::Array<long int>  MyCorrespondences;
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
      MyCorrespondences = pointvalue.Correspondence;
      scalars->InsertTuple1(idx, MyFA);
      tensors->SetTuple9(idx,MyTensor[0],MyTensor[1],MyTensor[2],MyTensor[3],MyTensor[4],MyTensor[5],MyTensor[6],MyTensor[7],MyTensor[8]);
      if (copyField.Correspondences)
      {
    	  correspondences->InsertTuple1(idx,MyCorrespondences[0]);
      }
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

    if (copyField.Correspondences)
    {
      polydata->GetPointData()->SetScalars(correspondences);
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
    correspondences->Delete();
  //}
  return polydata;
}

MeshType::Pointer vtk2itkMesh(vtkPolyData* polydata)
{
  MeshType::Pointer mesh = MeshType::New();
//  CopyVtkPolyData2ItkMesh(polydata, mesh);
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
        cellvalue.CaseName = caseName;
        popMesh->SetCell(citNew, copiedtrajectory);
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
