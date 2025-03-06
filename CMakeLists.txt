cmake_minimum_required(VERSION 3.15)

project(FlowThroughTensors LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)

set(SOURCES src/FlowThroughTensors.cpp)

add_library(FlowThroughTensors SHARED ${SOURCES})

# Include the current directory for headers
target_include_directories(FlowThroughTensors PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

if(OpenMP_CXX_FOUND)
    target_link_libraries(FlowThroughTensors PRIVATE OpenMP::OpenMP_CXX)
endif()

set_target_properties(FlowThroughTensors PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(WIN32)
    set_target_properties(FlowThroughTensors PROPERTIES OUTPUT_NAME "FlowThroughTensors")
elseif(APPLE)
    set_target_properties(FlowThroughTensors PROPERTIES OUTPUT_NAME "FlowThroughTensors_arm")
else()
    set_target_properties(FlowThroughTensors PROPERTIES OUTPUT_NAME "FlowThroughTensors")
endif()
