# Configure Project
cmake_minimum_required(VERSION 3.10)
project(Orbit_Simulator VERSION 1.0.0 LANGUAGES CXX)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set C++ standard (now you can use 17!)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include CTest module
include(CTest)
enable_testing()

# Find packages
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(jsoncpp REQUIRED)

# Generate executable
add_executable(Orbit_Simulator 
    main.cpp
    src/Simulator.cpp
)

# No need for BOOST_MATH_DISABLE_THREADS anymore!

# Link libraries
# Link libraries
target_link_libraries(Orbit_Simulator 
    PRIVATE 
    Eigen3::Eigen
    Boost::boost
    JsonCpp::JsonCpp  # <- Change this line
)

target_include_directories(Orbit_Simulator PRIVATE src)

# Move executable
add_custom_command(TARGET Orbit_Simulator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_FILE:Orbit_Simulator>
        ${CMAKE_SOURCE_DIR}/executable/Orbit_Simulator  # <- Remove .exe
    COMMAND ${CMAKE_COMMAND} -E remove 
        $<TARGET_FILE:Orbit_Simulator>    
)