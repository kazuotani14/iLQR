cmake_minimum_required(VERSION 3.5.0)
project (iLQR)

set (CMAKE_CXX_STANDARD 11)

###
# My stuff
###

# `cmake -DCMAKE_BUILD_TYPE=Debug ..` to compile with debug flags

# file(GLOB_RECURSE Foo_SOURCES "src/*.cpp")

set (Foo_SOURCES "src/run_ilqr.cpp"
                 "src/double_integrator.cpp"
                 "src/ilqr_core.cpp"
                 "src/forward_pass.cpp"
                 "src/derivatives.cpp"
                 "src/finite_diff.cpp"
             "src/boxqp.cpp"
         "src/eigen_helpers.cpp")
file(GLOB_RECURSE Foo_HEADERS "include/*.h")

set (Foo_INCLUDE_DIRS "")
foreach (_headerFile ${Foo_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND Foo_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES Foo_INCLUDE_DIRS)

include_directories(${Foo_INCLUDE_DIRS})
add_executable (run_iLQR ${Foo_SOURCES})


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
    add_executable(${exec_name} ${filepath})
    target_link_libraries(${exec_name} libgtest libgmock)
endmacro(add_gtest)

add_gtest(test_gtest "test/test_gtest.cpp")

add_executable(deriv "test/test_derivatives.cpp" "src/finite_diff.cpp")
target_link_libraries(deriv libgtest libgmock)
add_executable(boxqp "test/test_boxqp.cpp" "src/boxqp.cpp" "src/eigen_helpers.cpp")
target_link_libraries(boxqp libgtest libgmock)
add_executable(double_int "test/test_doubleintegrator.cpp" "src/double_integrator.cpp")
target_link_libraries(double_int libgtest libgmock)
