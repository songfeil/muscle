//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_GENERATE_MUSCLE_H
#define INTRODUCTION_GENERATE_MUSCLE_H

#include <Eigen/Core>
#include <vector>

void generate_muscle(const Eigen::MatrixXd & points,
                    std::vector<Eigen::MatrixXd> & VV,
                    std::vector<Eigen::MatrixXi> & FF);

#endif //INTRODUCTION_GENERATE_MUSCLE_H
