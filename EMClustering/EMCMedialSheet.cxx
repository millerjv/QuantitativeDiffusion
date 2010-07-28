#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
//// SLICER VERSION
#include "Common.h"
#include "EMClusteringIO.h"
#include "Registration.h"
#include "myMaths.h"
#include "EMCore.h"
#include "MeshOperations.h"
#include "Quantification.h"

#include "itkBinaryMask3DMeshSource.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkIdTypeArray.h>
#include "itkAutomaticTopologyMeshSource.h"
#include <itkQuadrilateralCell.h>


#include "EMCMedialSheetCLP.h"

ImageType::Pointer itkMeshToBinaryVolume(const MeshType* mesh, const ImageType* refVolume, int label)
{
	ImageType::Pointer vol = ImageType::New();
	ImageType::RegionType desiredRegion;
	desiredRegion.SetSize(  refVolume->GetLargestPossibleRegion().GetSize());
	desiredRegion.SetIndex( refVolume->GetLargestPossibleRegion().GetIndex());

	vol->SetRegions(desiredRegion);
	vol->Allocate();

	vol->SetSpacing(refVolume->GetSpacing());
	vol->SetOrigin(refVolume->GetOrigin());

	vol->FillBuffer(0);

	typedef MeshType::PointsContainerConstIterator pointsIt;
	pointsIt pointIt = mesh->GetPoints()->Begin();
	pointsIt end = mesh->GetPoints()->End();
	ImageType::IndexType ind;
    while (pointIt != end)
    {
    	if (vol->TransformPhysicalPointToIndex(pointIt.Value(), ind))
			{
			 vol->SetPixel(ind,label);
			}
    	else
    	{
    		std::cout << "point out of reference volume" <<std::endl;
    	}
    	pointIt++;

    }

	return vol;
}

ImageType* removeHoles(const ImageType* vol)
{
	typedef itk::VotingBinaryIterativeHoleFillingImageFilter< ImageType >  FilterType;
	FilterType::Pointer filter = FilterType::New();

	filter->SetInput(vol);

	const unsigned int radiusX = 2;
	const unsigned int radiusY = 2;
	const unsigned int radiusZ = 2;

    ImageType::SizeType indexRadius;

	indexRadius[0] = radiusX; // radius along x
	indexRadius[1] = radiusY; // radius along y
	indexRadius[2] = radiusZ; // radius along y

	filter->SetRadius( indexRadius );

	filter->SetBackgroundValue(   0 );
	filter->SetForegroundValue( 1 );
	filter->SetMajorityThreshold( 2 );

	filter->SetMaximumNumberOfIterations( 10 );
	filter->Update();

	return filter->GetOutput();
}

bool isEqual(std::vector<int> l1, std::vector<int> l2)
{
	bool flag=1;
	if (l1.size() != l2.size())
		flag = 0;
	else
	{
		for (unsigned int i=0; i<l1.size(); i++)
		{
			if (l1[i]!=l2[i])
			{
			  flag =0;
			  break;
			}
		}
	}
	return flag;
}


bool  hasNoTouchingFace(std::vector <std::vector<int> > linkLists, vtkIdList* ptsIds, int* facePIds)
{
	bool hasNoNei;
	std::vector<int> l1,l2,l3,l0;
	int id[4];
	id[0] = ptsIds->GetId(facePIds[0]);
	id[1] = ptsIds->GetId(facePIds[1]);
	id[2] = ptsIds->GetId(facePIds[2]);
	id[3] = ptsIds->GetId(facePIds[3]);
	// if the vertices are not unique, repeated link list should be avoided

	l0= linkLists.at(id[0]);
	if (id[1]!=id[0])
	{
		l1= linkLists.at(id[1]);
	}
	else
	{
		l1.push_back(-1);
	}
	if (id[2]!=id[0] && id[2]!=id[1])
	{
	    l2= linkLists.at(id[2]);
	}
	else
	{
		l2.push_back(-1);
	}
	if (id[3]!=id[0] && id[3]!=id[1] && id[3]!=id[2])
	{
		l3= linkLists.at(id[3]);
	}
	else
	{
		l3.push_back(-1);
	}
	if (isEqual(l0,l1) || isEqual(l0,l2) || isEqual(l0,l3) || isEqual(l1,l2) || isEqual(l1,l3) || isEqual(l2,l3))
	    hasNoNei = 0;
	else //unique
		hasNoNei = 1;

	return hasNoNei;
}

MeshType::PointType  ave(MeshType::PointType p1, MeshType::PointType p2)
{
	MeshType::PointType midPoint;
	midPoint[0]	= (p1[0]+p2[0])/2;
	midPoint[1]	= (p1[1]+p2[1])/2;
	midPoint[2]	= (p1[2]+p2[2])/2;
	return midPoint;
}



int main(int argc, char* argv[])
{
//	PARSE_ARGS;


	typedef itk::AutomaticTopologyMeshSource< MeshType > MeshSourceType;
	MeshSourceType::Pointer midSurface = MeshSourceType::New();

	//std::string filename ="/fs/corpus1/mahnaz/vtkMRMLFiniteElementMeshNode1.vtk";
	std::string filename ="/fs/corpus1/mahnaz/dtiCode/MORIAtlas_ROIs/IA-FEMesh-Autosave-20100728-1114/Model_4_4_VMesh-5.vtk";


    vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
    reader->SetFileName( filename.c_str() );
    std::cout<< "Reading "<<filename.c_str()<< "..."<<std::endl;
	reader->Update();
	vtkUnstructuredGrid* hexahedralGrid = reader->GetOutput();
	vtkIdType NCells = hexahedralGrid->GetNumberOfCells();
	hexahedralGrid->BuildLinks();
	vtkCellLinks* cellLinks = hexahedralGrid->GetCellLinks();

	vtkIdType NPoints = hexahedralGrid->GetNumberOfPoints();
	std::vector <std::vector<int> > linkLists;

	for (unsigned int pid=0; pid<NPoints; pid++ )
	{
	 vtkIdType *links = cellLinks->GetCells(pid);
	 int nlinks = cellLinks->GetNcells(pid);
	 std::vector<int> list;
	 for (int j=0; j<nlinks; j++)
	 {
		list.push_back(links[j]);
	 }
	 linkLists.push_back(list);
	}

	for (unsigned int c=0; c<NCells; c++)
	{
		vtkCell* aCell;
		aCell = hexahedralGrid->GetCell(c);
		vtkPoints* points = aCell->GetPoints();
		vtkIdList* ptsIds = aCell->GetPointIds();
		//face: {0,1,2,3}
		int face1PIds[4] = {0,1,2,3};
		bool face1 =  hasNoTouchingFace(linkLists, ptsIds, face1PIds);
        //face: {4,5,6,7}
		int face2PIds[4] = {4,5,6,7};
		bool face2 =  hasNoTouchingFace(linkLists, ptsIds, face2PIds);
		//face: {0,1,5,4}
		int face3PIds[4] = {0,1,5,4};
		bool face3 =  hasNoTouchingFace(linkLists, ptsIds, face3PIds);
		//face: {3,7,6,2}
		int face4PIds[4] = {3,7,6,2};
		bool face4 =  hasNoTouchingFace(linkLists, ptsIds, face4PIds);
		//face: {0,4,7,3}
		int face5PIds[4] = {0,4,7,3};
		bool face5 =  hasNoTouchingFace(linkLists, ptsIds, face5PIds);
		//face: {1,2,6,5}
		int face6PIds[4] = {1,2,6,5};
		bool face6 =  hasNoTouchingFace(linkLists, ptsIds, face6PIds);
		bool validMesh=1;
		std::cout << face1 << " " << face2 << " " << face3 << " " << face4 << " " <<face5 << " " << face6 <<std::endl;
		MeshType::PointType vp0, vp1, vp2, vp3;
		MeshType::PointType p0, p1, p2, p3;
		if (face3 && face4)  //{{0,1,5,4},{3,2,6,7}}
		{
			vp0 = ave(points->GetPoint(0), points->GetPoint(3));
			vp1 = ave(points->GetPoint(1), points->GetPoint(2));
			vp2 = ave(points->GetPoint(5), points->GetPoint(6));
			vp3 = ave(points->GetPoint(4), points->GetPoint(7));
		}
		else if (face1 && face2)  //{{0,1,2,3},{4,5,6,7}}
		{
			vp0 = ave(points->GetPoint(0), points->GetPoint(4));
			vp1 = ave(points->GetPoint(1), points->GetPoint(5));
			vp2 = ave(points->GetPoint(2), points->GetPoint(6));
			vp3 = ave(points->GetPoint(3), points->GetPoint(7));
		}

		else if (face5 && face6)  //{{0,3,7,4},{1,5,6,2}}
		{
			vp0 = ave(points->GetPoint(0), points->GetPoint(1));
			vp1 = ave(points->GetPoint(3), points->GetPoint(2));
			vp2 = ave(points->GetPoint(6), points->GetPoint(7));
			vp3 = ave(points->GetPoint(5), points->GetPoint(4));
		}
		else
		{
			validMesh = 0;
			std::cout << "Invalid input hexahedral mesh!" << std::endl;
		}
		//take care of itk vs. vtk difference
		if (validMesh)
		{
		p0[0] = -vp0[0];
		p0[1] = -vp0[1];
		p0[2] =  vp0[2];
		p1[0] = -vp1[0];
		p1[1] = -vp1[1];
		p1[2] =  vp1[2];
		p2[0] = -vp2[0];
		p2[1] = -vp2[1];
		p2[2] =  vp2[2];
		p3[0] = -vp3[0];
		p3[1] = -vp3[1];
		p3[2] =  vp3[2];
		//insert points and cells in the new mesh
		//std::cout << p0[0] << " " << p0[1] << " "<< p0[2] << std::endl;
		//std::cout << p1[0] << " " << p1[1] << " "<< p1[2] << std::endl;
		//std::cout << p2[0] << " " << p2[1] << " "<< p2[2] << std::endl;
		//std::cout << p3[0] << " " << p3[1] << " "<< p3[2] << std::endl;
		midSurface->AddQuadrilateral (midSurface->AddPoint(p1) , midSurface->AddPoint(p0),midSurface->AddPoint(p2), midSurface->AddPoint(p3));
		//int currentPoint =  midSurface->GetOutput()->GetNumberOfPoints()-1;
		//if (p0[2] >= 14.0212 || p1[2] >= 14.0212 || p2[2] >= 14.0212 || p3[2] >= 14.0212)
		if (c==275)
		{
			std::cout << c << std::endl;
			std::cout << ptsIds->GetId(0) << " " << ptsIds->GetId(1) << std::endl;
			std::cout << ptsIds->GetId(2) << " " << ptsIds->GetId(3) << std::endl;
			std::cout << ptsIds->GetId(4) << " " << ptsIds->GetId(5) << std::endl;
			std::cout << ptsIds->GetId(6) << " " << ptsIds->GetId(7) << std::endl;
			std::cout << face1 << " " << face2 << " " << face3 << " " << face4 << " " <<face5 << " " << face6 <<std::endl;
			std::cout << p0[0] << " " << p0[1] << " "<< p0[2] << std::endl;
			std::cout << p1[0] << " " << p1[1] << " "<< p1[2] << std::endl;
			std::cout << p2[0] << " " << p2[1] << " "<< p2[2] << std::endl;
			std::cout << p3[0] << " " << p3[1] << " "<< p3[2] << std::endl;
		}
		}
	}

	CopyFieldType copyfield ={0,0,0,0,0,0};
	WriteVTKfile(midSurface->GetOutput(), "/fs/corpus1/mahnaz/surface.vtp", copyfield);

    //midSurface->Print(std::cout);
    //midSurface->GetOutput()->Print(std::cout);
	reader->Delete();


	//Trajectories = ReadVTKfile(trajectoriesFilename.c_str());

    //Read FA map of the case
	//ImageType::Pointer caseFAVolume = ReadImageVolume(subjectFAFilename.c_str());

	//create a binary image volume
    //Trajectories_Volume = itkMeshToBinaryVolume(Trajectories, caseFAVolume, 1);

    //Trajectories_Volume = removeHoles(Trajectories_Volume);

    //Surface = extractSurface(Trajectories_Volume);

    return EXIT_SUCCESS;

}
