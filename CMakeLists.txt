cmake_minimum_required(VERSION 3.15)

project(neurodamus)

set(ND_CORE_DIR "" CACHE STRING "Additional Neurodamus mechanisms")
set(ND_ENABLE_CORENEURON ON CACHE BOOL "Enable CoreNEURON")
set(ND_MECHANISMS "" CACHE STRING "Which mechanisms to use")
set_property(CACHE ND_MECHANISMS PROPERTY STRINGS hippocampus neocortex thalamus)

option(ND_NCX_METABOLISM OFF)
option(ND_NCX_NGV OFF)
option(ND_NCX_PLASTICITY OFF)
option(ND_NCX_V5 OFF)

find_program(nrnivmodl nrnivmodl REQUIRED)

file(GLOB _input_core_hoc_files "${ND_CORE_DIR}/hoc/*.hoc")
file(GLOB _input_model_hoc_files "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/hoc/*.hoc")
set(input_hoc_files ${_input_core_hoc_files})
list(APPEND input_hoc_files ${_input_model_hoc_files})

set(input_mod_directories "${ND_CORE_DIR}/mod")

if("${ND_MECHANISMS}" STREQUAL "neocortex")
  if(ND_NCX_METABOLISM)
    list(APPEND input_mod_directories "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/mod/metabolism")
  endif()
  if(ND_NCX_NGV)
    list(APPEND input_mod_directories "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/mod/ngv")
  endif()
  if(ND_NCX_PLASTICITY)
    list(APPEND input_mod_directories "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/mod/v6/optimized")
  endif()
  list(APPEND input_mod_directories "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/mod/v6")
  if(ND_NCX_V5)
    list(APPEND input_mod_directories "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/mod/v5")
  endif()
else()
  list(APPEND input_mod_directories "${CMAKE_CURRENT_LIST_DIR}/${ND_MECHANISMS}/mod")
endif()

set(_input_mod_filestubs)

set(_archdir x86_64)
set(_build_mod_dir mods)

set(output_binaries ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/special)
set(output_libraries ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/libnrnmech.so)
set(output_cpp_files)
set(output_cpp_core_files)
set(output_mod_files)

if(ND_ENABLE_CORENEURON)
  set(_corenrn_arg "-coreneuron")
  list(APPEND output_binaries ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/special-core)
  list(APPEND output_libraries ${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/libcorenrnmech.so)
endif()

foreach(_input_dir IN LISTS input_mod_directories)
  file(GLOB _input_files "${_input_dir}/*.mod")
  foreach(_input_file IN LISTS _input_files)
    get_filename_component(_input_stub "${_input_file}" NAME_WLE)
    if("${_input_stub}" IN_LIST _input_mod_filestubs)
    else()
      set(_build_mod_file "${CMAKE_CURRENT_BINARY_DIR}/${_build_mod_dir}/${_input_stub}.mod")

      add_custom_command(COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${_input_file}" "${_build_mod_file}" OUTPUT "${_build_mod_file}")

      list(APPEND _input_mod_filestubs "${_input_stub}")
      list(APPEND output_cpp_files "${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/${_input_stub}.cpp")
      list(APPEND output_cpp_core_files "${CMAKE_CURRENT_BINARY_DIR}/${_archdir}/corenrn/mod2c/${_input_stub}.cpp")
      list(APPEND output_mod_files "${_build_mod_file}")
    endif()
  endforeach()
endforeach()

# At some point in the 1.14 series, HDF5 decided that the HDF5_*LIBRARIES variables should
# contain _targets_, not library names/paths.
set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
find_package(HDF5 REQUIRED)
find_package(sonata REQUIRED)

set(nrn_incflags "-isystem ${sonatareport_INCLUDE_DIR}")
set(nrn_loadflags "-Wl,-rpath,'${CMAKE_INSTALL_PREFIX}/lib' ${sonatareport_LIBRARY}")
foreach(dir IN LISTS HDF5_INCLUDE_DIRS)
  set(nrn_incflags "${nrn_incflags} -isystem '${dir}'")
endforeach()
foreach(lib IN LISTS HDF5_C_LIBRARIES)
  set(nrn_loadflags "${nrn_loadflags} ${lib}")
endforeach()

set(ENV{VERBOSE} "1")

add_custom_target(build ALL
  "${nrnivmodl}" ${_corenrn_arg} -incflags "${nrn_incflags}" -loadflags "${nrn_loadflags}" "${_build_mod_dir}"
  BYPRODUCTS ${output_binaries} ${output_libraries} ${output_cpp_files}
  DEPENDS ${output_mod_files}
)

INSTALL(PROGRAMS ${output_binaries} DESTINATION bin)
INSTALL(FILES ${output_libraries} DESTINATION lib)
INSTALL(FILES ${input_hoc_files} DESTINATION lib/hoc)
INSTALL(FILES ${output_mod_files} DESTINATION lib/mod)
INSTALL(FILES ${output_cpp_files} DESTINATION share/modcpp)

if(ND_ENABLE_CORENEURON)
  INSTALL(FILES ${output_cpp_core_files} DESTINATION share/modcpp_core)
endif()
