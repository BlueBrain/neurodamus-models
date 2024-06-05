cmake_minimum_required(VERSION 3.15)

project(neurodamus)

set(NEURODAMUS_CORE_DIR
    ""
    CACHE STRING "Base mechanisms needed for Neurodamus execution")
set(NEURODAMUS_ENABLE_CORENEURON
    ON
    CACHE BOOL "Enable CoreNEURON")
set(NEURODAMUS_MECHANISMS
    ""
    CACHE STRING "Which model's mechanisms to use")
set_property(CACHE NEURODAMUS_MECHANISMS PROPERTY STRINGS hippocampus neocortex
                                                  thalamus)

option(NEURODAMUS_NCX_METABOLISM OFF)
option(NEURODAMUS_NCX_NGV OFF)
option(NEURODAMUS_NCX_PLASTICITY OFF)
option(NEURODAMUS_NCX_V5 OFF)

# input_hoc_files will be installed via CMake's install(), which will basically
# just copy symbolic links. Make sure that we resolve any links ourselves!
set(input_hoc_files)
file(GLOB _input_core_hoc_files "${NEURODAMUS_CORE_DIR}/hoc/*.hoc")
file(GLOB _input_model_hoc_files
     "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/hoc/*.hoc")
foreach(_input_hoc_file IN LISTS _input_core_hoc_files _input_model_hoc_files)
  get_filename_component(_realpath "${_input_hoc_file}" REALPATH)
  list(APPEND input_hoc_files "${_realpath}")
endforeach()

# Mod directories in order of importance, files in later mod directories ignored
# if the same basename is present in an earlier directory
set(input_mod_directories "${NEURODAMUS_CORE_DIR}/mod")
if("${NEURODAMUS_MECHANISMS}" STREQUAL "neocortex")
  if(NEURODAMUS_NCX_METABOLISM)
    list(APPEND input_mod_directories
         "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/mod/metabolism")
  endif()
  if(NEURODAMUS_NCX_NGV)
    list(APPEND input_mod_directories
         "${CMAKE_CURRENT_LIST_DIR}/common/mod/ngv")
  endif()
  if(NEURODAMUS_NCX_PLASTICITY)
    list(APPEND input_mod_directories
         "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/mod/plasticity")
  endif()
  list(APPEND input_mod_directories
       "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/mod/v6")
  if(NEURODAMUS_NCX_V5)
    list(APPEND input_mod_directories
         "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/mod/v5")
  endif()
else()
  list(APPEND input_mod_directories
       "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/mod")
endif()

set(_input_mod_filestubs)
set(_build_mod_dir mod)

foreach(_input_dir IN LISTS input_mod_directories)
  message(STATUS "processing directory ${_input_dir}")
  file(GLOB _input_files "${_input_dir}/*.mod")
  foreach(_input_file IN LISTS _input_files)
    get_filename_component(_input_stub "${_input_file}" NAME_WLE)
    if("${_input_stub}" IN_LIST _input_mod_filestubs)
      message(STATUS "skipping ${_input_file}")
    else()
      message(STATUS "building ${_input_file}")
      list(APPEND mod_files ${_input_file})
      list(APPEND _input_mod_filestubs "${_input_stub}")
    endif()
  endforeach()
endforeach()

# NEURON does not install any CMake glue, we have to manually find nrnivmodl
# find_program(nrnivmodl nrnivmodl REQUIRED)
find_package(neuron REQUIRED)

# At some point in the 1.14 series, HDF5 decided that the HDF5_*LIBRARIES
# variables should contain _targets_, not library names/paths.
set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
find_package(HDF5 REQUIRED)
find_package(sonata REQUIRED)

if(NEURODAMUS_ENABLE_CORENEURON)
  create_nrnmech(CORENEURON SPECIAL INSTALL_MOD INSTALL_CPP MECHANISM_NAME neurodamus_${NEURODAMUS_MECHANISMS} MOD_FILES
                 ${mod_files})
  target_link_libraries(corenrnmech PRIVATE HDF5::HDF5 sonata::sonata_report)
  target_link_libraries(nrnmech PRIVATE HDF5::HDF5 sonata::sonata_report)
else()
  create_nrnmech(SPECIAL INSTALL_MOD INSTALL_CPP neurodamus_${NEURODAMUS_MECHANISMS} MOD_FILESMOD_FILES ${mod_files})
  target_link_libraries(nrnmech PRIVATE HDF5::HDF5 sonata::sonata_report)
endif()

install(FILES ${input_hoc_files}
        DESTINATION share/neurodamus_${NEURODAMUS_MECHANISMS}/hoc)

