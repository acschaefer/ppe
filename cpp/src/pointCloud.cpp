#include <pointCloud.h>

namespace ppe {

void pointCloud::resize(uint32_t width, uint32_t height) {
  this->width = width;
  this->height = height;
  Location.resize(width * height, 3);
  // TODO implement hasColor flag
  Color.resize(width * height, 3);
}

std::pair<uint32_t, uint32_t> pointCloud::subscript(uint32_t ip) {
  uint32_t row = ip / width;
  uint32_t col = ip % width;
  std::pair<uint32_t, uint32_t> sub(row, col);
  return sub;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud::createPCLcloud() {
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>());
  cloud->is_dense = false;
  // TODO check for color size
  cloud->width = width;
  cloud->height = height;
  cloud->resize(width*height);
  for (uint32_t i = 0; i < Location.rows(); i++) {
    cloud->points[i].x = Location(i, 0);
    cloud->points[i].y = Location(i, 1);
    cloud->points[i].z = Location(i, 2);
    cloud->points[i].r = Color(i, 0);
    cloud->points[i].g = Color(i, 1);
    cloud->points[i].b = Color(i, 2);
  }
  return cloud;
}

}
