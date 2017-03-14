cmake_minimum_required(VERSION 3.7.0)
project (iLQR)

include_directories(include)

file (GLOB SOURCES "src/*.cpp")

add_executable(run_iLQR src/main.cpp)
add_executable(test_dynamics src/loco_car.cpp)
