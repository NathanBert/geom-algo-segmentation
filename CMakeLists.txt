cmake_minimum_required(VERSION 3.5)
project(GeomAlgo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(Boost REQUIRED COMPONENTS headers)  # fixe warning

add_subdirectory(src)


target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE 
    "${CMAKE_SOURCE_DIR}/external/libigl/include"
    "${CMAKE_SOURCE_DIR}/external/libigl/external/glfw/include"  # si besoin
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
