cmake_minimum_required(VERSION 3.10)
project(TNWL VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(TNWL
        main.cpp
        input_graph.cpp
        tg.cpp)

target_link_libraries(TNWL PUBLIC voro++ libvoro++.a)
target_include_directories(TNWL PUBLIC
        "${PROJECT_BINARY_DIR}"
        "/usr/local/bin"
        "/usr/local/lib"
        "/usr/local/include"
        "/usr/local/include/voro++"
)
