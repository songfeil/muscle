//
// Created by jnyao on 11/7/18.
//

#ifndef INTRODUCTION_BEZIER_H
#define INTRODUCTION_BEZIER_H

#include "Eigen/Core"

void bezier(const Eigen::Vector3d & p0,
            const Eigen::Vector3d & p1,
            const Eigen::Vector3d & p2,
            const int n,
            Eigen::MatrixXd & B);

#endif //INTRODUCTION_BEZIER_H
