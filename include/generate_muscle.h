//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_GENERATE_MUSCLE_H
#define INTRODUCTION_GENERATE_MUSCLE_H

#include <Eigen/Core>
#include <vector>
#include <set>

//NEW ORGANIZATION

void generate_muscle(const Eigen::MatrixXd & points,
                     const int num_points,
                     Eigen::MatrixXd & Vm,
                     Eigen::MatrixXi & Fm);

void attach_muscle( const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    Eigen::MatrixXd & Vm,
                    Eigen::MatrixXi & Fm,
                    std::set<int> & attached_vids);

void attach_muscle_multiface( const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    Eigen::MatrixXd & Vm,
                    Eigen::MatrixXi & Fm);

void attach_tendon(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    const Eigen::MatrixXd & Vm,
                    const Eigen::MatrixXi & Fm,
                    std::vector<Eigen::MatrixXd> & VV,
                    std::vector<Eigen::MatrixXi> & FF);

#endif //INTRODUCTION_GENERATE_MUSCLE_H
