#ifndef POINTCLOUD_H_
#define POINTCLOUD_H_

#include <Eigen/Dense>

// TODO check for PCL
#include <pcl/io/pcd_io.h>
//#include <pcl/point_types.h>

namespace ppe {
class pointCloud {
public:
  uint32_t width;
  uint32_t height;
  Eigen::MatrixXd Location;
  Eigen::Matrix<uint8_t, -1, -1> Color;

  void resize(uint32_t width, uint32_t height);

  // TODO check for PCL
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr createPCLcloud();
};
}

#endif /* POINTCLOUD_H_ */
