cmake_minimum_required (VERSION 2.8.12)
project (flappie C)

option (BUILD_SHARED_LIB "Build a shared library" OFF)

set(CMAKE_CONFIGURATION_TYPES "Debug;Chaos;Release")
if(NOT CMAKE_BUILD_TYPE)
	message("Defaulting to release build")
	set(CMAKE_BUILD_TYPE Release)
endif()

include (InstallRequiredSystemLibraries)
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Scrappie is the technology demonstration platform for the Research algorithms' group.")
set (CPACK_PACKAGE_VENDOR "Oxford Nanopore Technologies")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.md")
#set (CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.md")
set (CPACK_PACKAGE_VERSION_MAJOR 1)
set (CPACK_PACKAGE_VERSION_MINOR 4)
set (CPACK_PACKAGE_VERSION_PATCH 0)

# Get the latest abbreviated commit hash of the working branch
execute_process(
	COMMAND git log -1 --format=%h
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_COMMIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-${GIT_COMMIT_HASH}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
set (CPACK_PACKAGE_NAME "ont-${PROJECT_NAME}")

set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Tim Massingham <tim.massingham@nanoporetech.com>")
set (CPACK_DEBIAN_PACKAGE_SECTION "base")
set (CPACK_DEBIAN_PACKAGE_DEPENDS "libopenblas-base, libhdf5-7, libcunit1")
set (CPACK_DEBIAN_BUILD_DEPENDS "libopenblas-base, libopenblas-dev, libhdf5-7, libhdf5-dev, cmake, libcunit1-dev")
set (CPACK_PACKAGING_INSTALL_PREFIX "/opt/flappie")

set (CPACK_GENERATOR "TGZ;DEB")
include (CPack)


configure_file (
    "${PROJECT_SOURCE_DIR}/src/version.h.in"
    "${PROJECT_BINARY_DIR}/include/version.h"
)


##
#   Set up include directories
##
include_directories ("${PROJECT_BINARY_DIR}/include")
if (OPENBLAS_ROOT)
	include_directories ("${OPENBLAS_ROOT}/include")
	link_directories ("${OPENBLAS_ROOT}/lib")
	set (BLAS "openblas")
else ()
	set (BLAS "blas")
endif ()
if (HDF5_ROOT)
	include_directories ("${HDF5_ROOT}/include")
	link_directories ("${HDF5_ROOT}/lib")
	set (CMAKE_REQUIRED_INCLUDES ${HDF5_ROOT}/include)
endif()


##
#   Set up what is to be built
##
add_library (flappie_objects OBJECT
	src/decode.c 
	src/layers.c 
	src/networks.c 
	src/nnfeatures.c 
	src/flappie_common.c 
	src/flappie_matrix.c 
        src/flappie_output.c
	src/util.c)
set_property(TARGET flappie_objects PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library (flappie_static STATIC $<TARGET_OBJECTS:flappie_objects>)
set_target_properties(flappie_static PROPERTIES OUTPUT_NAME flappie CLEAN_DIRECT_OUTPUT 1)
add_executable (flappie
	src/fast5_interface.c 
	src/flappie.c)

if (BUILD_SHARED_LIB)
	if (APPLE)
		message (SEND_ERROR "Building shared library on OSX not yet supported")
	endif (APPLE)
	add_library (flappie_shared SHARED $<TARGET_OBJECTS:flappie_objects>)
	set_target_properties(flappie_shared PROPERTIES OUTPUT_NAME flappie CLEAN_DIRECT_OUTPUT 1)
	install (TARGETS flappie_shared LIBRARY DESTINATION lib)
endif (BUILD_SHARED_LIB)

if (NOT DEFINED CHAOSMONKEY)
    set(CHAOSMONKEY 0.1)
endif ()


set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -Wall -Wunused-function -Wunused-value -Wunused-parameter -fstack-protector-all -fgnu89-inline -O3 -march=native -std=c99 -DUSE_SSE2 -D__USE_MISC -D_POSIX_SOURCE -DNDEBUG")
set (CMAKE_C_FLAGS_CHAOS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function -fstack-protector-all -fgnu89-inline -g -march=native -std=c99 -DUSE_SSE2 -D__USE_MISC -D_POSIX_SOURCE -DNDEBUG -DCHAOSMONKEY=${CHAOSMONKEY}")
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -Werror -Wall -Wno-cpp -DABORT_ON_NULL -Wno-unused-function -fstack-protector-all -fgnu89-inline -g -march=native -std=c99 -DUSE_SSE2 -D__USE_MISC -D_POSIX_SOURCE")


# Find right hdf5 file
include (CheckIncludeFile)
check_include_file ("hdf5.h" HDF5_STANDARD)
if (HDF5_STANDARD)
	set (HDF5 "hdf5")
else (HDF5_STANDARD)
	check_include_file ("hdf5/serial/hdf5.h" HDF5_SERIAL)
	if (HDF5_SERIAL)
		set (HDF5 "hdf5_serial")
		include_directories ("/usr/include/hdf5/serial")
	endif (HDF5_SERIAL)
endif (HDF5_STANDARD)

target_link_libraries (flappie flappie_static ${BLAS} ${HDF5} m)
if (APPLE)
	target_link_libraries (flappie argp)
endif (APPLE)

install (TARGETS flappie flappie_static RUNTIME DESTINATION bin ARCHIVE DESTINATION lib)


enable_testing()
add_executable(flappie_unittest 
	src/test/flappie_test_runner.c 
	src/test/flappie_util.c 
	src/test/test_flappie_convolution.c 
	src/test/test_flappie_elu.c 
	src/test/test_flappie_matrix.c 
	src/test/test_flappie_signal.c 
	src/test/test_flappie_util.c 
	src/test/test_skeleton.c 
	src/test/test_util.c)
target_include_directories(flappie_unittest PUBLIC "src/test" "src")
target_link_libraries(flappie_unittest flappie_static ${BLAS} ${HDF5} m cunit)

set (READSDIR ${PROJECT_SOURCE_DIR}/reads)
set (TESTREAD "MINICOL228_20161012_FNFAB42578_MN17976_mux_scan_HG_52221_ch271_read66_strand")
set (ENV{OPENBLAS_NUM_THREADS} 1)

add_test(NAME unittest WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/test/ COMMAND flappie_unittest)
add_test(test_call flappie ${READSDIR})
add_test(test_licence flappie --licence)
add_test(test_license flappie --license)
add_test(test_help flappie --help)
add_test(test_version flappie --version)

add_custom_target(test-verbose COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
