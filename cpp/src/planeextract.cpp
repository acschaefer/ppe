#include <iostream>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include "ptxread.h"

int main (int argc, char** argv)
{

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud = ple::ptxread("../../data/cottage_01_head.ptx");


  pcl::visualization::PCLVisualizer viewer ("PCL Viewer");

   // Define R,G,B colors for the point cloud
  //pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_color_handler(cloud, 255, 255, 255);
  // We add the point cloud to the viewer and pass the color handler
  //viewer.addPointCloud(cloud, cloud_color_handler, "cloud");
  viewer.addPointCloud(cloud, "cloud");

  viewer.addCoordinateSystem(1.0, "cloud", 0);
  viewer.setBackgroundColor(0.05, 0.05, 0.05, 0); // Setting background to a dark grey
  viewer.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");

  while (!viewer.wasStopped ()) { // Display the visualiser until 'q' key is pressed
    viewer.spinOnce ();
  }

  return 0;
}
