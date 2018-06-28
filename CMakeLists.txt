cmake_minimum_required(VERSION 2.8)
project(stella)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(cxx11)
include(cxx_config)
include(testing)

option(ENABLE_UNIT_TESTS "Enable unit testing" OFF)
option(AutodetectMPI "Automatically detect MPI" ON)

enable_language(C CXX Fortran)

find_package(PETSc REQUIRED)
if(PETSC_VERSION VERSION_LESS 3.7.0)
  message(FATAL_ERROR "ERROR: stella requires PETSc 3.7.0 or newer.")
endif()

set(req-lib ${PETSC_LIBRARIES} m)
set(req-inc ${PETSC_INCLUDES})

if(AutodetectMPI)
  find_package(MPI REQUIRED)
  list(APPEND req-lib ${MPI_C_LIBRARIES})
  list(APPEND req-inc ${MPI_C_INCLUDE_PATH})
endif()

include_directories(src)

if(ENABLE_UNIT_TESTS)
  add_subdirectory(test)
endif()

add_subdirectory(src)
add_library(stella ${stella_sources})
target_compile_definitions(stella PUBLIC PETSC_OPTIONS_TAKES_DATABASE)
target_link_libraries(stella PUBLIC ${req-lib})
target_include_directories(stella PUBLIC ${req-inc})

add_subdirectory(tools)

# Set the Doxygen list
set(list_for_doxygen_headers
    "src/stella_bc.h"
    "src/stella_boundary.h"
    "src/stella_classify.h"
    "src/stella_dirichlet.h"
    "src/stella_dmap.h"
    "src/stella_fd.h"
    "src/stella_gen.h"
    "src/stella_grid.h"
    "src/stella_io.h"
    "src/stella_level.h"
    "src/stella_mat.h"
    "src/stella_metric.h"
    "src/stella_neumann.h"
    "src/stella_operator.h"
    "src/stella_pc.h"
    "src/stella_rhs.h"
    "src/stella_signals.h"
    "src/stella_solver.h"
    "src/stella_state.h"
    "src/stella_stencil.h"
    "tools/single-block/base.h"
    "tools/single-block/block.h"
    "tools/single-block/boundary.h"
    "tools/single-block/grid.h"
    "tools/single-block/map.h"
    "tools/single-block/option.h"
    "tools/single-block/problem.h"
    "tools/single-block/solver.h"
    "tools/single-block/state.h"
    )
foreach(arg ${list_for_doxygen_headers})
    set(doxygen_headers "${doxygen_headers} ${CMAKE_CURRENT_SOURCE_DIR}/${arg}")
endforeach(arg ${list_for_doxygen_headers})

# build documentation (if enabled)
option(DOC "Create and install the HTML based API documentation (requires Doxygen)" ON)
if(DOC)

    find_package(Doxygen)

    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.bin)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION share/doc/stella OPTIONAL)
endif()
