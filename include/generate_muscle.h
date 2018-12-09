//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_GENERATE_MUSCLE_H
#define INTRODUCTION_GENERATE_MUSCLE_H

#include <Eigen/Core>
#include <vector>
#include <set>

void generate_muscle(const Eigen::MatrixXd & points,
                     const int num_points,
                     const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     const std::set<int> & selected_faces,
                     Eigen::MatrixXd & Vm,
                     Eigen::MatrixXi & Fm,
                     std::set<int> & fixed_vids);

void generate_muscle_multiface(const Eigen::MatrixXd & points,
                     const int num_points,
                     const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     const std::set<int> & selected_faces1,
                     const std::set<int> & selected_faces2,
                     std::vector<Eigen::MatrixXd> & VV,
                     std::vector<Eigen::MatrixXi> & FF);

#endif //INTRODUCTION_GENERATE_MUSCLE_H
