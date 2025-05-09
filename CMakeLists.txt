# Minimum required version of CMake
cmake_minimum_required(VERSION 3.15)

# Project name
project(wavimg LANGUAGES CXX)

if(WIN32)
    # !!!!! important !!!!!
    # For building on Windows, set this to the path of your vcpkg installation
    include(C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Where to find our headers
include_directories("${PROJECT_SOURCE_DIR}/include")

# Collect all .cpp files under src/
file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Add executable
add_executable(wavimg ${SOURCES})

find_package(cxxopts CONFIG REQUIRED)
target_link_libraries(wavimg PRIVATE cxxopts::cxxopts)

find_package(AudioFile CONFIG REQUIRED)
target_link_libraries(wavimg PRIVATE AudioFile)

find_package(lodepng CONFIG REQUIRED)
target_link_libraries(wavimg PRIVATE lodepng)

# Place final binary in <project_root>/bin
set_target_properties(wavimg PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
