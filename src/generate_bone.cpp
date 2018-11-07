#include "generate_bone.h"

void generate_bone(const Eigen::Vector3d &p0,
                   const Eigen::Vector3d &p1,
                   Eigen::MatrixXd &V,
                   Eigen::MatrixXi &F) {
  // Magic number
  int part = 5;

  Eigen::Vector3d vec = p1 - p0;
  Eigen::Vector3d u, v;
  u = Eigen::Vector3d(0, 1, vec(1) / vec(2)).normalized();
  v = u.cross(vec).normalized();



}