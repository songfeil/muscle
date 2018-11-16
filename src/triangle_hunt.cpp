//
// Created by 宋飞龙 on 2018-11-15.
//

#include "triangle_hunt.h"
#include <limits>

void triangle_hunt(Eigen::Matrix3d & P, Eigen::MatrixXd & V, Eigen::MatrixXd & F) {
  Eigen::RowVector3d triangleCenter = (P.row(0) + P.row(1) + P.row(2)) / 3.0;
  double dist = std::numeric_limits<double>::infinity();
  int fi = -1;

  for (int i = 0; i < F.rows(); i++) {
    Eigen::RowVector3i f = F.row(i);
    Eigen::RowVector3d p0 = V.row(f(0));
    Eigen::RowVector3d p1 = V.row(f(1));
    Eigen::RowVector3d p2 = V.row(f(2));

    Eigen::RowVector3d pc = (p0 + p1 + p2) / 3.0;

    double cdist = (pc - triangleCenter).norm();

    if (cdist < dist) {
      fi = i;
    }
  }

  // Change the location
  Eigen::RowVector3i f = F.row(fi);
  V.row(f(0)) = P.row(0);
  V.row(f(1)) = P.row(1);
  V.row(f(2)) = P.row(2);
}