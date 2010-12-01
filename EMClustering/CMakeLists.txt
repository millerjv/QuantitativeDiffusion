PROJECT(Clustering)

set(CMAKE_BUILD_TYPE Debug CACHE STRING "Configuration type being built (needed to run tests).")
enable_testing()

cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# --------------------------------------------------------------------------
# Find Slicer3

if(NOT Slicer3_SOURCE_DIR)
  find_package(Slicer3 REQUIRED)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
endif(NOT Slicer3_SOURCE_DIR)

# Include dirs

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${Slicer3_Libs_INCLUDE_DIRS}
  ${Slicer3_Base_INCLUDE_DIRS}
)

set (CLP EMClustering)

# --------------------------------------------------------------------------
# Build and install the library

set(lib_name ClusteringLib)
add_library(${lib_name}
    EMClusteringIO.cxx AffineRegistration.cxx myMaths.cxx EMCore.cxx MeshOperations.cxx Quantification.cxx QuadEdge/itkQuadEdge.cxx)
  
slicer3_set_modules_output_path(${lib_name})
target_link_libraries(${lib_name}
  ${Slicer3_Libs_LIBRARIES}
  ${Slicer3_Base_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
)

set_target_properties (${lib_name} PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
slicer3_set_plugins_output_path(${lib_name})
slicer3_install_plugins(${lib_name})

GenerateCLP(${CLP}.cxx ${CLP}.xml)
add_executable(${CLP} ${CLP}.cxx)
target_link_libraries(${CLP} ${lib_name})
slicer3_set_plugins_output_path(${CLP})

slicer3_install_plugins(${CLP})

#-----------testing  
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



