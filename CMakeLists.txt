# Configure Project
cmake_minimum_required(VERSION 3.10)
project(Orbit_Simulator VERSION 1.0.0 LANGUAGES CXX)

set(BUILD_TESTING OFF CACHE BOOL "Disable testing" FORCE)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find packages - Use Conan's target names
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(jsoncpp REQUIRED)

# Build the nrlmsise-00 library
add_subdirectory(dependencies/nrlmsise-00)

# Generate executable
add_executable(Orbit_Simulator 
    main.cpp
    src/Simulator.cpp
)

# Link libraries - Use lowercase target names for Conan packages
target_link_libraries(Orbit_Simulator 
    PRIVATE 
    Eigen3::Eigen
    boost::boost
    JsonCpp::JsonCpp
    nrlmsise00  
)

target_include_directories(Orbit_Simulator PRIVATE src)

# Copy executable
if(WIN32)
    add_custom_command(TARGET Orbit_Simulator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
            $<TARGET_FILE:Orbit_Simulator>
            ${CMAKE_SOURCE_DIR}/executable/Orbit_Simulator.exe
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:nrlmsise00>
            ${CMAKE_SOURCE_DIR}/executable/
        COMMAND ${CMAKE_COMMAND} -E remove 
            $<TARGET_FILE:Orbit_Simulator>
    )
else()
    # For Linux
    add_custom_command(TARGET Orbit_Simulator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
            $<TARGET_FILE:Orbit_Simulator>
            ${CMAKE_SOURCE_DIR}/executable/Orbit_Simulator
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:nrlmsise00>
            ${CMAKE_SOURCE_DIR}/executable/
        COMMAND ${CMAKE_COMMAND} -E remove 
            $<TARGET_FILE:Orbit_Simulator>
    )
endif()