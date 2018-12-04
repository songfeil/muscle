//
// Created by 宋飞龙 on 2018-12-03.
//

#ifndef INTRODUCTION_MOVE_PATCH_H
#define INTRODUCTION_MOVE_PATCH_H

#include <Eigen/Core>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>

Eigen::Vector3d calc_center(Eigen::MatrixXd & V) {
  Eigen::Vector3d sum(0, 0, 0);
  for (int i = 0; i < V.rows(); i++) {
    sum += V.row(i).transpose();
  }
  sum = sum / V.rows();
  return sum;
}

void translate_V(Eigen::MatrixXd & V, const Eigen::Vector3d & t) {
  for (int i = 0; i < V.rows(); i++) {
    V.row(i) += t.transpose();
  }
}

void origin_scaling(Eigen::MatrixXd & V, double s) {
  Eigen::Vector3d center = calc_center(V);
  translate_V(V, -center);
  V = V * s;
  translate_V(V, center);
}

void move_patch(
    Eigen::MatrixXd & musV,
    Eigen::MatrixXi & musF,
    Eigen::MatrixXd & bonV,
    Eigen::MatrixXi & bonF
    ) {
  Eigen::Vector3d musCenter = calc_center(musV);
  Eigen::Vector3d bonCenter = calc_center(bonV);

  Eigen::Vector3d direction = (bonCenter - musCenter).normalized();

  igl::Hit hit;
  bool result = igl::ray_mesh_intersect(musCenter, direction, bonV, bonF, hit);
  if (!result) exit(4);

  origin_scaling(musV, 1.2);
  translate_V(musV, hit.t * direction);
}

#endif //INTRODUCTION_MOVE_PATCH_H