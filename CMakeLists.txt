cmake_minimum_required(VERSION 3.15)

project(neurodamus)

set(NEURODAMUS_CORE_DIR "" CACHE STRING "Base mechanisms needed for Neurodamus execution")
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

# NEURON does not install any CMake glue, we have to manually find nrnivmodl
find_program(nrnivmodl nrnivmodl REQUIRED)

# Extract the nrnivmodl output directory
file(STRINGS "${nrnivmodl}" _modsubdir REGEX "MODSUBDIR=.*")
string(REPLACE "MODSUBDIR=" "" _archdir "${_modsubdir}")
message(STATUS "nrnivmodl output directory determined to be \"${_archdir}\"")

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
         "${CMAKE_CURRENT_LIST_DIR}/${NEURODAMUS_MECHANISMS}/mod/v6/optimized")
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

set(output_binaries ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/special)
set(output_libraries
    ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/libnrnmech${CMAKE_SHARED_LIBRARY_SUFFIX}
)
set(output_cpp_files)
set(output_cpp_core_files)
set(output_mod_files)

if(NEURODAMUS_ENABLE_CORENEURON)
  set(_corenrn_arg "-coreneuron")
  list(APPEND output_binaries
       ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/special-core)
  list(
    APPEND
    output_libraries
    ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/libcorenrnmech${CMAKE_SHARED_LIBRARY_SUFFIX}
  )
endif()

foreach(_input_dir IN LISTS input_mod_directories)
  file(GLOB _input_files "${_input_dir}/*.mod")
  foreach(_input_file IN LISTS _input_files)
    get_filename_component(_input_stub "${_input_file}" NAME_WLE)
    if("${_input_stub}" IN_LIST _input_mod_filestubs)
      message(STATUS "skipping ${_input_file}")
    else()
      message(STATUS "building ${_input_file}")
      set(_build_mod_file
          "${CMAKE_CURRENT_BINARY_DIR}/${_build_mod_dir}/${_input_stub}.mod")

      add_custom_command(
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${_input_file}"
                "${_build_mod_file}" OUTPUT "${_build_mod_file}")

      list(APPEND _input_mod_filestubs "${_input_stub}")
      list(APPEND output_cpp_files
           "${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/${_input_stub}.cpp")
      list(
        APPEND
        output_cpp_core_files
        "${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/corenrn/mod2c/${_input_stub}.cpp"
      )
      list(APPEND output_mod_files "${_build_mod_file}")
    endif()
  endforeach()
endforeach()

# At some point in the 1.14 series, HDF5 decided that the HDF5_*LIBRARIES
# variables should contain _targets_, not library names/paths.
set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)
find_package(sonata REQUIRED)

set(nrn_incflags "-isystem ${sonatareport_INCLUDE_DIR}")
get_filename_component(_report_dir "${sonatareport_LIBRARY}" DIRECTORY)
set(nrn_loadflags "-Wl,-rpath,'${_report_dir}' '${sonatareport_LIBRARY}'")
foreach(dir IN LISTS HDF5_INCLUDE_DIRS MPI_C_INCLUDE_DIRS)
  set(nrn_incflags "${nrn_incflags} -isystem '${dir}'")
endforeach()
foreach(lib IN LISTS HDF5_C_LIBRARIES MPI_C_LIBRARIES)
  get_filename_component(_lib_dir "${lib}" DIRECTORY)
  set(nrn_loadflags "${nrn_loadflags} -Wl,-rpath,'${_lib_dir}' '${lib}'")
endforeach()

configure_file(build_neurodamus.sh.in build_neurodamus.sh @ONLY)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(nrn_loadflags "-Wl,-rpath,'${CMAKE_INSTALL_PREFIX}/lib' ${nrn_loadflags}")

add_custom_target(
  build ALL
  "${nrnivmodl}"
  ${_corenrn_arg}
  -incflags
  "${nrn_incflags}"
  -loadflags
  "${nrn_loadflags}"
  "${_build_mod_dir}"
  BYPRODUCTS ${output_binaries} ${output_libraries} ${output_cpp_files}
  DEPENDS ${output_mod_files})

install(PROGRAMS ${output_binaries}
                 ${CMAKE_CURRENT_BINARY_DIR}/build_neurodamus.sh
        DESTINATION bin)
install(FILES ${output_libraries} DESTINATION lib)
install(FILES ${input_hoc_files}
        DESTINATION share/neurodamus_${NEURODAMUS_MECHANISMS}/hoc)
install(FILES ${output_mod_files}
        DESTINATION share/neurodamus_${NEURODAMUS_MECHANISMS}/mod)
install(FILES ${output_cpp_files}
        DESTINATION share/neurodamus_${NEURODAMUS_MECHANISMS}/cpp)

if(NEURODAMUS_ENABLE_CORENEURON)
  install(FILES ${output_cpp_core_files}
          DESTINATION share/neurodamus_${NEURODAMUS_MECHANISMS}/cpp_core)
endif()
