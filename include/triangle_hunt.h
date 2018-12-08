//
// Created by 宋飞龙 on 2018-11-15.
//

#ifndef INTRODUCTION_TRIANGLE_HUNT_H
#define INTRODUCTION_TRIANGLE_HUNT_H

#include <Eigen/Core>

//P is a 3 * 3 matrix containing a mesh triangle, where every row is position of vertex.

void triangle_hunt(Eigen::Matrix3d & P, Eigen::MatrixXd & V, Eigen::MatrixXi & F);
int triangle_hunts(Eigen::Matrix3d & P, Eigen::MatrixXd V, Eigen::MatrixXi F);
void triangle_hunt_lst(Eigen::RowVector3d & P, const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::MatrixXd & lst_dist);

#endif //INTRODUCTION_TRIANGLE_HUNT_H
