#ifndef PTXREAD_H_
#define PTXREAD_H_

#include <Eigen/Dense>

#include <pointCloud.h>

namespace ple {
pointCloud* ptxread(std::string filename);
}

#endif /* PTXREAD_H_ */
