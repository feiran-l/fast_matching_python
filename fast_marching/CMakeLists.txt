cmake_minimum_required(VERSION 3.14)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE "Release" )
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC")


if(APPLE)
    message("------- The current platform is macos! ---------")
    set(CMAKE_CXX_COMPILER clang++)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
elseif(UNIX AND NOT APPLE)
    message("------- The current platform is linux! ---------")
    set(CMAKE_CXX_FLAGS "-fopenmp -Wall -Wno-unused-result")
endif()

project(fast_marching LANGUAGES CXX)


find_package(OpenMP REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE CPP_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(fast_marching_lib STATIC ${CPP_SOURCES})
target_link_libraries(fast_marching_lib ${OpenMP_CXX_LIBRARIES})
target_link_libraries(fast_marching_lib Eigen3::Eigen)
set_target_properties(fast_marching_lib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})


add_subdirectory(pybind11)
pybind11_add_module(fast_marching fm_module.cpp)
target_link_libraries(fast_marching PRIVATE fast_marching_lib)
set_target_properties(fast_marching PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# cmake -DPYTHON_EXECUTABLE:FILEPATH=python3 ..






