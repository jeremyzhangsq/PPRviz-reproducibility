cmake_minimum_required(VERSION 3.10)
project(approx_dnppr)

set(CMAKE_CXX_STANDARD 11)
#SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -DNDEBUG")

#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Boost REQUIRED system serialization)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# download + cmake + make install + include directories: http://eigen.tuxfamily.org/index.php?title=Main_Page
find_package(Eigen3)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}")
file(GLOB SOURCES "*.cpp")
add_executable(approx_dnppr ${SOURCES})

find_package(OpenMP REQUIRED)
if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP_TARGET INTERFACE)
    add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
    target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
    target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)

