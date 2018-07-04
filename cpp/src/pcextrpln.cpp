#include <pcextrpln.h>

using namespace Eigen;

namespace ppe {

pointCloud* pcextrpln(pointCloud* pc) {

  // Define the different initial plane configurations.
  std::vector<Eigen::Matrix3i> c;

  // Compute the radii of all rays.
  Eigen::MatrixXd l(pc->Location); // TODO redundant copy
  Eigen::VectorXd r(l.rowwise().norm());

  // Compute the normalized ray direction vectors.
  Eigen::MatrixXd v(l.rowwise().normalized());

  //std::cout << pc->Location << std::endl;
  //std::cout << r << std::endl;
  //std::cout << v << std::endl;

  //Compute the plane parameters and the errors corresponding to creating
  // planes at different locations and in different configurations.
  for (uint32_t ip = 0; ip < l.rows(); ip++) {
    // Determine the subscript index of the center point of the plane.
    std::pair<uint32_t, uint32_t> oc = pc->subscript(ip);

  }

  return pc;
}

}
