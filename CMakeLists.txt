cmake_minimum_required(VERSION 3.10)
project(ProjectileTrajectory)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(EXECUTABLE_NAME main)

file(GLOB SOURCES "src/*.cpp")  # Adjust path/pattern as needed

# Capture both .hpp and .h headers (in include/ and src/ here)
file(GLOB HEADERS
    "include/*.hpp"
    "include/*.h"
    "src/*.hpp"
    "src/*.h"
)

# Create the executable
add_executable(${EXECUTABLE_NAME} ${SOURCES} ${HEADERS})

# Include directories (if you have headers in other locations)
target_include_directories(${EXECUTABLE_NAME} PRIVATE include/)  # Adjust as needed