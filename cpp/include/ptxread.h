#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace ple {
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr ptxread(std::string filename);
}
