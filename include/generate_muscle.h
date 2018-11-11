//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_GENERATE_MUSCLE_H
#define INTRODUCTION_GENERATE_MUSCLE_H

#include <Eigen/Core>

void generate_bone(const Eigen::MatrixXd & B,
                   const Eigen::MatrixXd & N,
                   Eigen::MatrixXd & V,
                   Eigen::MatrixXi & F);

#endif //INTRODUCTION_GENERATE_MUSCLE_H
