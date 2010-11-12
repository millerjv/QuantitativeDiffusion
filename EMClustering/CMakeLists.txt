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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(ClusteringLib EMClusteringIO.cxx AffineRegistration.cxx myMaths.cxx EMCore.cxx MeshOperations.cxx Quantification.cxx QuadEdge/itkQuadEdge.cxx) # SurfaceOperations.cxx)
set_target_properties (ClusteringLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
target_link_libraries(ClusteringLib ITKCommon ITKIO vtkIO vtkCommon vtkFiltering vtkGraphics vtkRendering ITKAlgorithms)

set (CLP_LIST EMClustering ) #EMClustering2010 EMMedialSheetExtract transform_vtk vtp2vtk)

foreach(CLP ${CLP_LIST})
set ( ${CLP}_SOURCE ${CLP}.cxx)
GenerateCLP(${CLP}.cxx ${CLP}.xml)

add_executable(${CLP} ${CLP}.cxx)
target_link_libraries( ${CLP} ClusteringLib ITKCommon ITKIO vtkIO vtkCommon vtkFiltering vtkGraphics vtkRendering ITKAlgorithms)
  
slicer3_set_plugins_output_path(${CLP})
endforeach(CLP)

set(Slicer3_EXE ${Slicer3_DIR}/Slicer3)

if (WIN32)
  set( PluginsPath ${Clustering_BINARY_DIR}/lib/Slicer3/Plugins/${CMAKE_BUILD_TYPE}/ )
else (WIN32)
  set( PluginsPath ${Clustering_BINARY_DIR}/lib/Slicer3/Plugins/)
endif (WIN32)

add_test(EMClusteringTest1 ${PluginsPath}/EMClustering
${Clustering_BINARY_DIR}/Testing/Temporary/OutputClusters.vtp 
--trajectories ${Clustering_SOURCE_DIR}/Data/trajectories_test2.vtp --initialCenters ${Clustering_SOURCE_DIR}/Data/initCenters_test2.vtp
--outdirectory ${Clustering_BINARY_DIR}/Testing/Temporary)


# configure the resources for the atlas
set(ATLAS_PARENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB ATLASFILES RELATIVE ${ATLAS_PARENT_DIR} "Atlas/?*.*")
foreach(file ${ATLASFILES} )
    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_LIB_DIR}/${PROJECT_NAME}/${file} COPYONLY)
endforeach(file)

install(FILES ${ATLASFILES} DESTINATION ${Slicer3_INSTALL_PLUGINS_LIB_DIR}/${PROJECT_NAME}/Atlas)

slicer3_install_plugins(EMClustering)
