cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(breaked VERSION 1.0)

add_subdirectory(src/)

target_include_directories(breaked PUBLIC "${breaked_SOURCE_DIR}/include")

# compile libraries
add_subdirectory(third-party/spdlog)
add_subdirectory(third-party/csv-parser)
add_subdirectory(third-party/argparse)
