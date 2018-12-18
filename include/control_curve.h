//
// Created by jnyao on 11/7/18.
//

#ifndef INTRODUCTION_BEZIER_H
#define INTRODUCTION_BEZIER_H

#include "Eigen/Core"

// Generate a Catmull Rom Spline curve given the control points p and
// output nsample points along the curve in C as well as their
// respective normals in N
void CatmullRomChain(const Eigen::MatrixXd & p,
            const int nsample,
            Eigen::MatrixXd & C,
            Eigen::MatrixXd & N);

#endif //INTRODUCTION_BEZIER_H
