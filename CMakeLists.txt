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

add_library(ekf_core STATIC
    src/tools.cpp
    src/dataset_parser.cpp
    src/FusionEKF.cpp
    src/kalman_filter.cpp
)

target_include_directories(ekf_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(ekf_core PUBLIC
    Eigen3::Eigen
)

target_compile_features(ekf_core PUBLIC cxx_std_17)

if(BUILD_TESTING)
    find_package(Catch2 3 CONFIG REQUIRED)

    add_executable(ekf_tests
        tests/test_tools.cpp
        tests/test_dataset_parser.cpp
        tests/test_FusionEKF.cpp
        tests/test_kalman_filter.cpp
    )

    target_link_libraries(ekf_tests PRIVATE
        ekf_core
        Catch2::Catch2WithMain
    )

    include(Catch)
    catch_discover_tests(ekf_tests)
endif()

# Apply compiler warnings only to maintained project targets.
set(ekf_warning_targets ekf_core)

if(BUILD_TESTING)
    list(APPEND ekf_warning_targets ekf_tests)
endif()

foreach(target IN LISTS ekf_warning_targets)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
        )
    endif()
endforeach()
