cmake_minimum_required(VERSION 3.10)
project(ProjectileTrajectory)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(COMMON_SOURCES
    src/derivative_functions.cpp
    src/simulation.cpp
)

add_executable(main src/main.cpp ${COMMON_SOURCES})
add_executable(compare_integrators src/compare_integrators.cpp ${COMMON_SOURCES})

# Include directories
target_include_directories(main PRIVATE include)
target_include_directories(compare_integrators PRIVATE include)
