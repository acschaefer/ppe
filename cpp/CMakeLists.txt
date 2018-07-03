cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(planeextract)

find_package(PCL 1.7 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_subdirectory(src)
