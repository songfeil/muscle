//
// Created by jnyao on 11/7/18.
//

#ifndef INTRODUCTION_BEZIER_H
#define INTRODUCTION_BEZIER_H

#include "Eigen/Core"

void bezier(const Eigen::MatrixXd & p,
//            const Eigen::Vector3d & p1,
//            const Eigen::Vector3d & p2,
            const int ndegree,
            const int nsample,
            Eigen::MatrixXd & B,
            Eigen::MatrixXd & N);

void CatmullRomChain(const Eigen::MatrixXd & p,
//            const Eigen::Vector3d & p1,
//            const Eigen::Vector3d & p2,
            const int nsample,
            Eigen::MatrixXd & C,
            Eigen::MatrixXd & N);

#endif //INTRODUCTION_BEZIER_H
