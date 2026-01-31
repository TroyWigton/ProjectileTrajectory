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
add_executable(test src/test.cpp ${COMMON_SOURCES})
add_executable(compare_k_over_m src/compare_k_over_m.cpp ${COMMON_SOURCES})
add_executable(test_variable_drag src/test_variable_drag.cpp ${COMMON_SOURCES})
add_executable(example_3d_motion src/example_3d_motion.cpp)
add_executable(test_golf_ball src/test_golf_ball.cpp ${COMMON_SOURCES})
add_executable(test_adaptive src/test_adaptive.cpp ${COMMON_SOURCES})

# Include directories
target_include_directories(main PRIVATE include)
target_include_directories(compare_integrators PRIVATE include)
target_include_directories(test PRIVATE include)
target_include_directories(compare_k_over_m PRIVATE include)
target_include_directories(test_variable_drag PRIVATE include)
target_include_directories(example_3d_motion PRIVATE include)
target_include_directories(test_golf_ball PRIVATE include)
target_include_directories(test_adaptive PRIVATE include)
