cmake_minimum_required(VERSION 3.15)

project(FlowThroughTensors LANGUAGES CXX)

# C++14 is essential for generic lambdas used for sorting in the source code
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

add_executable(${PROJECT_NAME} src/NetFlow/FlowThroughTensors.cpp)
