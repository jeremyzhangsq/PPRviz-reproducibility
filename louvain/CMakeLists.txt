cmake_minimum_required(VERSION 3.10)
project(louvainplus)

set(CMAKE_CXX_STANDARD 11)


set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")



find_package(Boost)

IF (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
endif()

find_package(OpenMP REQUIRED)
if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP_TARGET INTERFACE)
    add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
    target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
    target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()

include_directories(. ../approx_dnppr/)

add_executable(louvainplus
        community.cpp
        community.h
        graph.cpp
        graph.h
        graph_binary.cpp
        graph_binary.h
        main_community.cpp)

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)