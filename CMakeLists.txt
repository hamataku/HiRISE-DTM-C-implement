cmake_minimum_required(VERSION 3.10)

project(hirise_cpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")

find_package(GDAL REQUIRED)
if (GDAL_FOUND)
    include_directories(${GDAL_INCLUDE_DIR})
    set(GDAL_LIBS ${GDAL_LIBRARY})
endif()

add_executable(main main.cpp)

target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(main PRIVATE ${GDAL_LIBS})
