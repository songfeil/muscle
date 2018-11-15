//
// Created by 宋飞龙 on 2018-11-11.
//

#include "deform.h"
#include <Eigen/Geometry>
#include <math.h>
#include "pick_constrain_point.h"
#include "bezier.h"
#include <igl/arap.h>
#include "generate_bone.h"

void deform(
    const Eigen::Vector3d &p0,
    const Eigen::Vector3d &p1,
    const Eigen::Vector3d &p2,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F
) {
    std::vector<int> b;
    pick_constrain_point(V, b);

//    // Move
//    Eigen::Vector3d vec = p1 - p0;
//    double len = vec.norm();
//    double orilen = 5.0;
//    double s = len / orilen;
//    Eigen::MatrixXd scale = Eigen::MatrixXd(3, 3);
//    scale << s, 0, 0,
//             0, s, 0,
//             0, 0, s;
//
//    Eigen::Vector3d y_norm = Eigen::Vector3d(0, 1, 0);
//    double costheta = vec.normalized().dot(y_norm);
//    double theta = acos(costheta);
//    Eigen::Matrix3d rot = Eigen::Matrix3d(3, 3);
//    rot << cos(theta), -sin(theta), 0,
//           sin(theta), cos(theta), 0,
//           0, 0, 1;
//
//    for (int i = 0; i < V.rows(); i++) {
//        // Scale, Rotate and translate.
//        V.row(i) = V.row(i) * scale * rot + p0.transpose();
//    }

    generate_bone(p0, p1, V, F);

    // Deformation
//    Eigen::VectorXi bb = Eigen::VectorXi::Map(b.data(), b.size());;
    Eigen::VectorXi bb = Eigen::VectorXi(5);
    bb << 0, 4, 8, 12, 16;
    Eigen::MatrixXd Bc, N;
    bezier(p0, p2, p1, 4, Bc, N);
    igl::ARAPData data;
    igl::arap_precomputation(V, F, 3, bb, data);
    igl::arap_solve(Bc, data, V);

}