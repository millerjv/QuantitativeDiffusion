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
/*
MeshType::Pointer extractSurface(const ImageType* vol)
{
	MeshType::Pointer surface = MeshType::New();
	typedef itk::BinaryMask3DMeshSource< ImageType, MeshType >   MeshSourceType;
    MeshSourceType::Pointer meshSource = MeshSourceType::New();
    meshSource->SetObjectValue( 1 );
    meshSource->SetInput(vol);
    meshSource->Update();
    surface = meshSource->GetOutput();
	return surface;
}
*/
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

int main(int argc, char* argv[])
{
	PARSE_ARGS;

	MeshType::Pointer    Trajectories, Surface;
	ImageType::Pointer Trajectories_Volume;

	//Trajectories = ReadVTKfile(trajectoriesFilename.c_str());

    //Read FA map of the case
	//ImageType::Pointer caseFAVolume = ReadImageVolume(subjectFAFilename.c_str());

	//create a binary image volume
    //Trajectories_Volume = itkMeshToBinaryVolume(Trajectories, caseFAVolume, 1);

    //Trajectories_Volume = removeHoles(Trajectories_Volume);

    //Surface = extractSurface(Trajectories_Volume);

    return EXIT_SUCCESS;

}
