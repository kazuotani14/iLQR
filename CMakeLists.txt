cmake_minimum_required(VERSION 3.5.0)
project (iLQR)

set (CMAKE_CXX_STANDARD 11)

###
# My stuff
###

# `cmake -DCMAKE_BUILD_TYPE=Debug ..` to compile with debug flags

# file(GLOB_RECURSE Foo_SOURCES "src/*.cpp")

set (ILQR_SOURCES "src/double_integrator.cpp"
                  "src/ilqr_core.cpp"
                  "src/forward_pass.cpp"
                  "src/derivatives.cpp"
                  "src/finite_diff.cpp"
                  "src/boxqp.cpp"
                  "src/eigen_helpers.cpp"
                  "src/backward_pass.cpp")

file(GLOB_RECURSE ILQR_HEADERS "include/*.h")

set (ILQR_INCLUDE_DIRS "")
foreach (_headerFile ${ILQR_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND ILQR_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES ILQR_INCLUDE_DIRS)

include_directories(${ILQR_INCLUDE_DIRS})
add_executable (run_iLQR "src/run_ilqr.cpp" ${ILQR_SOURCES})


###
# Use ExternalProject module to fetch GoogleTest and GoogleMock from internet and compile as third-party dependencies.
# From: http://kaizou.org/2014/11/gtest-cmake/
###

# We need thread support
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Download and install GoogleTest
ExternalProject_Add(
    gtest
    URL https://github.com/google/googletest/archive/master.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # Disable install step
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(libgtest PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# Create a libgmock target to be used as a dependency by test programs
add_library(libgmock IMPORTED STATIC GLOBAL)
add_dependencies(libgmock gtest)

# Set libgmock properties
set_target_properties(libgmock PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
include_directories("${source_dir}/googletest/include"
                    "${source_dir}/googlemock/include")

###
# My tests
###

# TODO add optional arguments
macro(add_gtest exec_name filepath)
    set(OtherFiles ${ARGN})
    add_executable(${exec_name} ${filepath} ${OtherFiles})
    target_link_libraries(${exec_name} libgtest libgmock)
endmacro(add_gtest)

add_gtest(gtest_test "test/test_gtest.cpp")

add_gtest(deriv_test "test/test_finite_diff.cpp" "src/finite_diff.cpp")
add_gtest(boxqp_test "test/test_boxqp.cpp" "src/boxqp.cpp" "src/eigen_helpers.cpp")
add_gtest(double_int_test "test/test_doubleintegrator.cpp" "src/double_integrator.cpp")
add_gtest(ilqr_deriv_test "test/test_ilqr_derivatives.cpp" ${ILQR_SOURCES})
add_gtest(ilqr_forward_test "test/test_ilqr_forward_pass.cpp" ${ILQR_SOURCES})
