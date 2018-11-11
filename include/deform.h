//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_DEFORM_H
#define INTRODUCTION_DEFORM_H

#include <Eigen/Core>

void deform(
      const Eigen::Vector3d &p0,
      const Eigen::Vector3d &p1,
      const Eigen::Vector3d &p2,
      Eigen::MatrixXd & V,
      Eigen::MatrixXi & F
    );

#endif //INTRODUCTION_DEFORM_H
