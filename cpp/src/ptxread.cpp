#include "ptxread.h"
#include <iostream>
#include <fstream>

namespace ple {
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr ptxread(std::string filename) {

    std::cout << "Reading PCL file: " << filename << std::endl;

    // TODO check valid file

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>());
    cloud->is_dense = false;

    std::ifstream file(filename.c_str());
    std::string line;
    long long il = 0;
    double intensity;
    unsigned int r, g, b;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      if (il == 0) {
	iss >> cloud->width;
      } else if (il == 1) {
	iss >> cloud->height;
	cloud->points.resize (cloud->width * cloud->height);
      } else if (il >= 10) {
	iss >> cloud->points[il-10].x >> cloud->points[il-10].y >> cloud->points[il-10].z;
	iss >> intensity >> r >> g >> b;
	//iss >> cloud->points[il-10].r >> cloud->points[il-10].g >> cloud->points[il-10].b;
	cloud->points[il-10].r = (uint8_t) r;
	cloud->points[il-10].g = (uint8_t) g;
	cloud->points[il-10].b = (uint8_t) b;
      }
      il++;
    }


    /*
    // Fill in the cloud data
    cloud->width    = 500;
    cloud->height   = 1;

    for (size_t i = 0; i < cloud->points.size (); ++i) {
      cloud->points[i].x = 1024 * rand () / (RAND_MAX + 1.0f);
      cloud->points[i].y = 1024 * rand () / (RAND_MAX + 1.0f);
      cloud->points[i].z = 1024 * rand () / (RAND_MAX + 1.0f);
    }
    */

    return cloud;
  }
}
