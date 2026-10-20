cmake_minimum_required(VERSION 3.22)

project(
    extended_kalman_filter
    VERSION 0.1.0
    DESCRIPTION "Lidar and Radar sensor fusion using Extended Kalman Filter"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

find_package(Eigen3 CONFIG REQUIRED)

if(BUILD_TESTING)
    find_package(Catch2 3 CONFIG REQUIRED)
endif()

