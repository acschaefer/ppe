#include <iostream>
#include <fstream>

#include <ptxread.h>

namespace ppe {
pointCloud* ptxread(std::string filename) {

  std::cout << "INFO: Reading PCL file: " << filename << std::endl;

  // TODO check valid file

  pointCloud* pc = new pointCloud();

  std::ifstream file(filename.c_str());
  std::string line;
  uint32_t il = 0;
  double x, y, z, intensity;
  unsigned int r, g, b;
  uint32_t width, height;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    if (il == 0) {
      iss >> width;
    } else if (il == 1) {
      iss >> height;
      // TODO check availability of intensity / color
      pc->resize(width, height);
    } else if (il >= 10 && il-10 >= width*height) {
      std::cout << "WARNING: Found more points than width*height. Stopping file reading." << std::endl;
      break;
    } else if (il >= 10) {
      iss >> x >> y >> z;
      iss >> intensity >> r >> g >> b;
      pc->Location(il-10, 0) = x;
      pc->Location(il-10, 1) = y;
      pc->Location(il-10, 2) = z;
      pc->Color(il-10, 0) = r;
      pc->Color(il-10, 1) = g;
      pc->Color(il-10, 2) = b;
    }
    il++;
  }
  return pc;
}
}
