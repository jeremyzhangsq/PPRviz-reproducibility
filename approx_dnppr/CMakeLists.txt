cmake_minimum_required(VERSION 3.10)
project(approx_dnppr)

if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
    endif()

endif()

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


target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)

