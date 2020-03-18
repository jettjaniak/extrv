cmake_minimum_required(VERSION 3.14)
project(extrv)
set(CMAKE_CXX_STANDARD 17)
include_directories(headers)

add_subdirectory(pybind11)

pybind11_add_module(extrv_engine wrappers.cpp
        sources/helpers.cpp
        sources/Ligand.cpp
        sources/SimulationState.cpp
        sources/velocities.cpp
        sources/forces.cpp
        sources/interpolated.cpp
        sources/Parameters.cpp
        sources/History.cpp
        sources/AbstractBondType.cpp)

add_executable(extrv main.cpp
        sources/helpers.cpp
        sources/Ligand.cpp
        sources/SimulationState.cpp
        sources/velocities.cpp
        sources/forces.cpp
        sources/interpolated.cpp
        sources/Parameters.cpp
        sources/History.cpp
        sources/AbstractBondType.cpp)


