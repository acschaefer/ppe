cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

include_directories(${PROJECT_SRC_DIR}/include)

add_executable (matrix_test matrix_test.cpp)
target_link_libraries (matrix_test ${PCL_LIBRARIES})

add_executable (pcd_write pcd_write.cpp)
target_link_libraries (pcd_write ${PCL_LIBRARIES})

add_library(ppe SHARED
  inc_test.cpp
  pcdraw.cpp
  pointCloud.cpp
  ptxread.cpp
)

add_executable (planeextract planeextract.cpp)
target_link_libraries (planeextract ${PCL_LIBRARIES} ppe)

add_executable (eigen_test eigen_test.cpp)
