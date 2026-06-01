cmake_minimum_required(VERSION 3.10)
project(ProjectileTrajectory)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(COMMON_SOURCES
    src/derivative_functions.cpp
    src/simulation.cpp
)

add_executable(optimal_angle_search src/optimal_angle_search.cpp ${COMMON_SOURCES})
add_executable(compare_integrators src/compare_integrators.cpp ${COMMON_SOURCES})
add_executable(test src/test.cpp ${COMMON_SOURCES})
add_executable(compare_k_over_m src/compare_k_over_m.cpp ${COMMON_SOURCES})
add_executable(mach_drag_study src/mach_drag_study.cpp ${COMMON_SOURCES})
add_executable(example_3d_motion src/example_3d_motion.cpp)
add_executable(golf_ball_drag_study src/golf_ball_drag_study.cpp ${COMMON_SOURCES})
add_executable(demo_rk45_adaptive src/demo_rk45_adaptive.cpp ${COMMON_SOURCES})
add_executable(convergence_study src/convergence_study.cpp ${COMMON_SOURCES})

# Include directories
target_include_directories(optimal_angle_search PRIVATE include)
target_include_directories(compare_integrators PRIVATE include)
target_include_directories(test PRIVATE include)
target_include_directories(compare_k_over_m PRIVATE include)
target_include_directories(mach_drag_study PRIVATE include)
target_include_directories(example_3d_motion PRIVATE include)
target_include_directories(golf_ball_drag_study PRIVATE include)
target_include_directories(demo_rk45_adaptive PRIVATE include)
target_include_directories(convergence_study PRIVATE include)
