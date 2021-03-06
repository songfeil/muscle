//
// Created by jnyao on 12/5/18.
//

#ifndef INTRODUCTION_PCA_ELIPSE_H
#define INTRODUCTION_PCA_ELIPSE_H

#include <Eigen/Core>
#include <math.h>
#include <igl/sort.h>
#include <igl/slice.h>
#include <Eigen/Eigenvalues>
#include <iostream>

// Struct that combines the parameters defining an ellipse
struct ellipse_param
{
    double long_axis;
    double short_axis;
    Eigen::Vector3d long_dir;
    Eigen::Vector3d short_dir;
    Eigen::MatrixXd center;
    Eigen::MatrixXd normal;
};

// Given a patch mesh specified by Vpatch and Fpatch, fit an ellipse
// using PCA by setting the first two principal directions as long and short axis
// and the last principal direction as normal
struct ellipse_param PCA_param(Eigen::MatrixXd Vpatch, Eigen::MatrixXi Fpatch) {
    struct ellipse_param result;
    Eigen::MatrixXd center = Vpatch.colwise().mean();
    result.center = center;
    Eigen::MatrixXd P(Vpatch.rows(), 3);
    for (int i = 0; i < Vpatch.rows(); i++) {
        P.row(i) = Vpatch.row(i) - center;
    }

    // Eigen decomposition on P^T * P
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(P.transpose() * P);
    Eigen::MatrixXd eigenval(es.eigenvalues().rows(), 1);
    eigenval.col(0) = es.eigenvalues();
    Eigen::MatrixXd Y, IX;
    igl::sort(eigenval, 1, false, Y, IX); // descending

    // Sort the eigenvalue eigenvector pairs
    Eigen::Matrix3d sortedVec;
    sortedVec.col(0) = es.eigenvectors().col(IX(0));
    sortedVec.col(1) = es.eigenvectors().col(IX(1));
    sortedVec.col(2) = es.eigenvectors().col(IX(2));

    // Long axis
    result.long_dir = sortedVec.col(0);
    result.long_axis = (P * result.long_dir).cwiseAbs().maxCoeff();

    // Short axis
    result.short_dir = sortedVec.col(1);
    result.short_axis = (P * result.short_dir).cwiseAbs().maxCoeff();

    // Normal
    result.normal = sortedVec.col(2);

    return result;
}

#endif //INTRODUCTION_PCA_ELIPSE_H
