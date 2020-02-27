cmake_minimum_required(VERSION 3.14)
project(extrv)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(pybind11)

include_directories(headers)
add_executable(extrv
        main.cpp
        sources/helpers.cpp
        sources/Ligand.cpp
        sources/SimulationState.cpp
        sources/velocities.cpp
        sources/forces.cpp
        sources/interpolated.cpp
        sources/Settings.cpp)

pybind11_add_module(engine engine.cpp)
set_target_properties(engine PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})