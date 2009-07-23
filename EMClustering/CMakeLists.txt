cmake_minimum_required(VERSION 2.4)
# This project is designed to be built outside the Insight source tree.
PROJECT(Clustering)

# Slicer3
find_package(Slicer3 REQUIRED)
include(${Slicer3_USE_FILE})


# Find GenerateCLP.
FIND_PACKAGE(GenerateCLP REQUIRED)
IF (GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ENDIF (GenerateCLP_FOUND)



# Default install prefix
slicer3_set_default_install_prefix_for_external_projects()


# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# Find VTK.
FIND_PACKAGE(VTK)
IF(VTK_FOUND)
	INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)


IF (GenerateCLP_FOUND)
  GenerateCLP(EMClustering.cxx EMClustering.xml)
  ADD_EXECUTABLE(EMClustering EMClustering.cxx)
  slicer3_set_plugins_output_path(EMClustering)
ENDIF (GenerateCLP_FOUND)

TARGET_LINK_LIBRARIES( EMClustering ITKCommon ITKIO VTKIO VTKCommon VTKFiltering VTKGraphics VTKRendering)


