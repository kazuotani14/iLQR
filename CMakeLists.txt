cmake_minimum_required(VERSION 3.7.0)
project (iLQR)

include_directories(include)

file (GLOB SOURCES "src/*.cpp")

# find_package(Eigen REQUIRED)

add_executable(run_iLQR src/loco_car_main.cpp
												src/loco_car_dynamics.cpp
												src/loco_car_costs.cpp
												src/loco_car_main.cpp
												src/iLQR_dynamics_derivatives.cpp
												src/iLQR_core.cpp
												src/iLQR_forward_pass.cpp)

#add_executable(run_loco src/loco_car_dynamics.cpp
#												src/loco_car_costs.cpp
#												src/loco_car_main.cpp
#												src/iLQR_dynamics_derivatives.cpp)
