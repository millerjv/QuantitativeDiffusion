project(QuantitativeDiffusionTools)

cmake_minimum_required(VERSION 2.4)

find_package(Slicer3 REQUIRED NO_DEFAULT_PATH)

if (Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

  include_directories( ${SlicerBaseCLI_BINARY_DIR} ${SlicerBaseCLI_SOURCE_DIR} )

  link_directories(${ModuleDescriptionParser_BINARY_DIR})
  link_directories(${SlicerBaseCLI_BINARY_DIR})

  add_subdirectory(EMClustering)

else (Slicer3_FOUND)

message(WARNING "Cannot build with a Slicer3 build tree or a Slicer3 installation. Please set Slicer3_DIR. If using a Slicer3 build, Slicer3_DIR needs to be set to the root of the build tree of Slicer3.  If using an installation of Slicer3, Slicer3_DIR needs to be set to <SlicerInstallation>/lib/Slicer3.")

endif (Slicer3_FOUND)

