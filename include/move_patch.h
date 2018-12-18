//
// Created by 宋飞龙 on 2018-12-03.
//

#ifndef INTRODUCTION_MOVE_PATCH_H
#define INTRODUCTION_MOVE_PATCH_H

#include <Eigen/Core>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <math.h>

// Calculate the center point for given vertices
Eigen::Vector3d calc_center(Eigen::MatrixXd & V) {
    Eigen::Vector3d sum(0, 0, 0);
    for (int i = 0; i < V.rows(); i++) {
        sum += V.row(i).transpose();
    }
    sum = sum / V.rows();
    return sum;
}

// Translate vertices V by vector t
void translate_V(Eigen::MatrixXd & V, const Eigen::Vector3d & t) {
    for (int i = 0; i < V.rows(); i++) {
        V.row(i) += t.transpose();
    }
}

// Uniform scaling for vertices V by s
void origin_scaling(Eigen::MatrixXd & V, double s) {
    Eigen::Vector3d center = calc_center(V);
    translate_V(V, -center);
    V = V * s;
    translate_V(V, center);
}

// Search for the given rad angle is between which two indices of angle lists
int between_vertices(double query_rad, const Eigen::VectorXd & rads) {
    // All rads between 2 * PI - 0
    int size = rads.size();
    for (int i = 0; i < size; i++) {
        int j = (i + 1) % size;
        int ti = (i + 1) > size ? 1 : 0;

        if (query_rad <= rads(i) && query_rad > rads(j) - ti * 2 * M_PI) {
            return i;
        }
    }

//    std::cout << "return -1\n";

    return size - 1;
}

// Move patch using translate method
void translate_move_patch(
        Eigen::MatrixXd & musV,
        Eigen::MatrixXi & musF,
        Eigen::MatrixXd & bonV,
        Eigen::MatrixXi & bonF
) {
    Eigen::Vector3d musCenter = calc_center(musV);
    Eigen::Vector3d bonCenter = calc_center(bonV);

    Eigen::Vector3d direction = (bonCenter - musCenter).normalized();

    igl::Hit hit;
    bool result = igl::ray_mesh_intersect(musCenter, direction, bonV, bonF, hit);
    if (!result) exit(4);

    origin_scaling(musV, 1.2);
    translate_V(musV, hit.t * direction);
}

// Move patch by parameterization
void map_move_patch(
        Eigen::MatrixXd & musV,
        Eigen::MatrixXi & musF,
        Eigen::MatrixXd & bonV,
        Eigen::MatrixXi & bonF,
        Eigen::MatrixXd & musVNew
) {
    Eigen::VectorXi bonBND_orig, bonBND;
    igl::boundary_loop(bonF, bonBND_orig);
    bonBND.resize(bonBND_orig.size());

    // Flip the boundary point order
    int start_idx = 0;
    double ndist = std::numeric_limits<double>::infinity();
    for (int i = 0; i < bonBND.size(); i++) {
        int idx = (start_idx + i) % bonBND.size();
        bonBND(i) = bonBND_orig(idx);
    }
    Eigen::VectorXi musBND_orig, musBND;
    igl::boundary_loop(musF, musBND_orig);
    musBND.resize(musBND_orig.size());


    // Flip the boundary point order
    for (int i = 0; i < musBND_orig.rows(); i++) {
        musBND(i) = musBND_orig(musBND_orig.rows()-1-i);
    }
    musBND_orig = musBND;

    Eigen::RowVector3d bonCenter = bonV.colwise().mean();
    Eigen::RowVector3d musCenter = musV.colwise().mean();
    // Pick the proper starting point in the loop
    start_idx = 0;
    ndist = std::numeric_limits<double>::infinity();
    Eigen::RowVector3d target = musV.row(musBND(0)) - musCenter;
    for (int i = 0; i < musBND.size(); i++) {
        Eigen::Vector3d q = musV.row(musBND_orig(i));
        double dot_product = - target(0) * q(0) + target(1) * q(1) + target(2) * q(2);
        if ((dot_product / (target.norm() * q.norm())) < ndist) {
            ndist = dot_product / (target.norm() * q.norm());
            start_idx = i;
        }
    }

    for (int i = 0; i < musBND.size(); i++) {
        int idx = (start_idx + i) % musBND.size();
        musBND(i) = musBND_orig(idx);
    }

    Eigen::MatrixXd bon_uv;
    igl::map_vertices_to_circle(bonV, bonBND, bon_uv);
    Eigen::MatrixXd bonU, bonBND_U;
    igl::harmonic(bonV,bonF,bonBND,bon_uv,1,bonU);

    Eigen::MatrixXd mus_uv;
    igl::map_vertices_to_circle(musV, musBND, mus_uv);
    Eigen::MatrixXd musU, musBND_U;
    igl::harmonic(musV,musF,musBND,mus_uv,1,musU);

//    std::cout << "bon_boundary" << std::endl;
//    for (int i = 0; i < bonBND.size(); i++) {
//        std::cout<<bonV.row(bonBND(i))<<std::endl;
//    }
//    std::cout << "mus_boundary" << std::endl;
//    for (int i = 0; i < musBND.size(); i++) {
//        std::cout<<musV.row(musBND(i)) <<std::endl;
//    }

    bonBND_U.resize(bonBND.size(), 2);
    for (int i = 0; i < bonBND.size(); i++) {
        bonBND_U.row(i) = bonU.row(bonBND(i));
    }

    musBND_U.resize(musBND.size(), 2);
    for (int i = 0; i < musBND.size(); i++) {
        musBND_U.row(i) = musU.row(musBND(i));
    }

    Eigen::VectorXd mus_bound_rad(musBND.size());
    for (int i = 0; i < musBND.size(); i++) {
        double rad = atan2(musBND_U(i,0), musBND_U(i, 1)) + M_PI;
        rad = rad < 0 ? rad + 2 * M_PI : rad;
        rad = rad + 0.5 * M_PI;
        rad = rad > 2 * M_PI ? rad - 2 * M_PI : rad;
        mus_bound_rad(i) = rad;
    }

    Eigen::VectorXd bon_bound_rad(bonBND.size());
    for (int i = 0; i < bonBND.size(); i++) {
        double rad = atan2(bonBND_U(i,0), bonBND_U(i, 1)) + M_PI;
        rad = rad < 0 ? rad + 2 * M_PI : rad;
        rad = rad + 0.5 * M_PI;
        rad = rad > 2 * M_PI ? rad - 2 * M_PI : rad;
        bon_bound_rad(i) = rad;

    }

//    std::cout << "mus_bound_rad" <<std::endl;
//    std::cout << mus_bound_rad <<std::endl;


//    std::cout << "info loaded" << std::endl;

    musVNew = musV;
    for (int i = 0; i < mus_bound_rad.size(); i++) {
        double q = mus_bound_rad(i);

        int idx = between_vertices(q, bon_bound_rad);
        int idy = (idx + 1) % bon_bound_rad.size();

//        std::cout << "query" << q << " between " << bon_bound_rad(idx) << "and" << bon_bound_rad(idy) << std::endl;

        double brady = idy == 0 ? bon_bound_rad(idy) - 2 * M_PI : bon_bound_rad(idy);

        double interval = bon_bound_rad(idx) - brady;
        double query_interval = bon_bound_rad(idx) - q;
        double percentage = query_interval / interval;

//        std::cout << "percentage " << percentage <<std::endl;

        Eigen::RowVector3d new_pos = (1 - percentage) * bonV.row(idx) + (percentage) * bonV.row(idy);
        musVNew.row(musBND(i)) = new_pos;
    }
//    std::cout << "cal finished" << std::endl;

}

#endif //INTRODUCTION_MOVE_PATCH_H