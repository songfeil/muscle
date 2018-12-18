//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_GENERATE_MUSCLE_H
#define INTRODUCTION_GENERATE_MUSCLE_H

#include <Eigen/Core>
#include <vector>
#include <set>

// Generate muscle mesh given control points
// Vm, Fm are the output muscle mesh
void generate_muscle(const Eigen::MatrixXd & points,
                     const int num_points,
                     Eigen::MatrixXd & Vm,
                     Eigen::MatrixXi & Fm);

// Attach muscle mush towards selected faces on bone
// selected_faces record the face index selected on the bone mesh F
// V, F define bone mesh
// Vm, Fm are input muscle mesh and will be modified
void attach_muscle( const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    Eigen::MatrixXd & Vm,
                    Eigen::MatrixXi & Fm,
                    std::set<int> & attached_vids);

// Attach muscle mush towards selected faces on bone
// selected_faces record the face index selected on the bone mesh F
// V, F define bone mesh
// Vm, Fm are input muscle mesh and will be modified
void attach_muscle_multiface( const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    Eigen::MatrixXd & Vm,
                    Eigen::MatrixXi & Fm,
                    std::set<int> & attached_vids);

// Connect muscle mush with selected faces on bone
// selected_faces record the face index selected on the bone mesh F
// V, F define bone mesh
// Vm, Fm are input muscle mesh
// VV, FF are vectors recording all the meshes to be displayed in the window
// and the output tendon mesh are pushed into VV and FF
void attach_tendon(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    const Eigen::MatrixXd & Vm,
                    const Eigen::MatrixXi & Fm,
                    std::vector<Eigen::MatrixXd> & VV,
                    std::vector<Eigen::MatrixXi> & FF,
                    std::set<int> & attached_vids);

#endif //INTRODUCTION_GENERATE_MUSCLE_H
