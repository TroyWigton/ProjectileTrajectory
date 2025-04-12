cmake_minimum_required(VERSION 3.10)
project(ProjectileTrajectory)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(EXECUTABLE_NAME main)  # Change this to your desired executable name

file(GLOB SOURCES "src/*.cpp")  # Adjust path/pattern as needed
file(GLOB HEADERS "include/*.h")  # Adjust path/pattern as needed

# Create the executable
add_executable(${EXECUTABLE_NAME} ${SOURCES} ${HEADERS})

# Include directories (if you have headers in other locations)
target_include_directories(${EXECUTABLE_NAME} PRIVATE include/)  # Adjust as needed