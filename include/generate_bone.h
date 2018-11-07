#ifndef INTRODUCTION_GENERATE_BONE_H
#define INTRODUCTION_GENERATE_BONE_H

#include <Eigen/Core>

//Generate an easy bone mesh given two points in space.
//Input:
//  p0  point 0
//  p1  point 1
//Output
//  V
//  F
void generate_bone(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                   Eigen::MatrixXd & V, Eigen::MatrixXi & F);

#endif //INTRODUCTION_GENERATE_BONE_H