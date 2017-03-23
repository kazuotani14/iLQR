cmake_minimum_required(VERSION 3.7.0)
project (iLQR)


# find_package(Eigen REQUIRED)

file(GLOB_RECURSE Foo_SOURCES "src/*.cpp")

file(GLOB_RECURSE Foo_HEADERS "include/*.h")

set (Foo_INCLUDE_DIRS "")
foreach (_headerFile ${Foo_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND Foo_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES Foo_INCLUDE_DIRS)

include_directories(${Foo_INCLUDE_DIRS})
add_executable (run_iLQR ${Foo_SOURCES})
