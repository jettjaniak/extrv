cmake_minimum_required(VERSION 3.14)
project(extrv)

set(CMAKE_CXX_STANDARD 14)

add_executable(extrv main.cpp headers/constants.h sources/helpers.cpp headers/helpers.h headers/SimulationState.h sources/Ligand.cpp sources/SimulationState.cpp sources/velocities.cpp headers/velocities.h sources/forces.cpp headers/forces.h headers/interpolated.h sources/interpolated.cpp sources/Settings.cpp)