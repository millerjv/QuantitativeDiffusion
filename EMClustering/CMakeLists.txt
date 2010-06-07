PROJECT(Clustering)

set(CMAKE_BUILD_TYPE Debug CACHE STRING "Configuration type being built (needed to run tests).")

cmake_minimum_required(VERSION 2.4)
enable_testing()

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

set (CLP EMClustering)
set ( ${CLP}_SOURCE ${CLP}.cxx)


GenerateCLP(${CLP}.cxx ${CLP}.xml)
ADD_EXECUTABLE(${CLP} ${CLP}.cxx)
TARGET_LINK_LIBRARIES( ${CLP} ITKCommon ITKIO vtkIO vtkCommon vtkFiltering vtkGraphics vtkRendering ITKAlgorithms)
  
add_library(${CLP}Lib SHARED ${${CLP}_SOURCE})
set_target_properties (${CLP}Lib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
target_link_libraries (${CLP}Lib ITKCommon ITKIO vtkIO vtkCommon vtkFiltering vtkGraphics vtkRendering ITKAlgorithms)
  
slicer3_set_plugins_output_path(${CLP})
 

set(Slicer3_EXE ${Slicer3_DIR}/Slicer3)

if (WIN32)
  set( PluginsPath ${Clustering_BINARY_DIR}/lib/Slicer3/Plugins/${CMAKE_BUILD_TYPE}/ )
else (WIN32)
  set( PluginsPath ${Clustering_BINARY_DIR}/lib/Slicer3/Plugins/)
endif (WIN32)

add_test(EMClusteringTest1 ${Slicer3_EXE} --launch ${PluginsPath}/EMClustering
${Clustering_SOURCE_DIR}/Data/trajectories_test2.vtp 
${Clustering_BINARY_DIR}/Testing/Temporary/OutputClusters.vtp 
--initialCenters ${Clustering_SOURCE_DIR}/Data/initCenters_test2.vtp
--outdirectory ${Clustering_BINARY_DIR}/Testing/Temporary)


#slicer3_install_plugins(EMClustering)

# configure the resources for the atlas
set(ATLAS_PARENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB ATLASFILES RELATIVE ${ATLAS_PARENT_DIR} "Atlas/?*.*")
foreach(file ${ATLASFILES} )
    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_LIB_DIR}/${PROJECT_NAME}/${file} COPYONLY)
endforeach(file)

install(FILES ${ATLASFILES} DESTINATION ${Slicer3_INSTALL_PLUGINS_LIB_DIR}/${PROJECT_NAME}/Atlas)
