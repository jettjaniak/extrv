cmake_minimum_required(VERSION 3.14)
project(extrv)

set(CMAKE_CXX_STANDARD 14)

add_executable(extrv main.cpp constants.h helpers.cpp helpers.h SimulationState.h Ligand.cpp SimulationState.cpp velocities.cpp velocities.h forces.cpp forces.h interpolated.h interpolated.cpp)