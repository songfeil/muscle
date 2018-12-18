//
// Created by jnyao on 11/7/18.
//

#include "control_curve.h"
#include <math.h>
#include <iostream>
#include <igl/cat.h>

void CatmullRomSpline(Eigen::MatrixXd P0, Eigen::MatrixXd P1, Eigen::MatrixXd P2, Eigen::MatrixXd P3, Eigen::MatrixXd & C, Eigen::MatrixXd & N,  int nPoints=10) {
    double alpha = 0.5;
    auto tj = [&](double ti, Eigen::MatrixXd Pi, Eigen::MatrixXd Pj) {
        double xi = Pi(0, 0);
        double yi = Pi(0, 1);
        double xj = Pj(0, 0);
        double yj = Pj(0, 1);
        return pow(pow(((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi)), 0.5), alpha) + ti;
    };

    double t0 = 0;
    double t1 = tj(t0, P0, P1);
    double t2 = tj(t1, P1, P2);
    double t3 = tj(t2, P2, P3);
    C.resize(nPoints, 3);
    N.resize(nPoints, 3);

    auto P1P2 = [&] (double t, Eigen::MatrixXd & Cc, int rowi) {
        Eigen::MatrixXd A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1;
        Eigen::MatrixXd A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2;
        Eigen::MatrixXd A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3;

        Eigen::MatrixXd B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2;
        Eigen::MatrixXd B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3;

        Cc.row(rowi) = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2;
    };

    for (int i = 0; i < nPoints; i++) {
        double t = t1 + (t2 - t1) / (nPoints) * i + (t2 - t1) / (2*nPoints);

        P1P2(t, C, i);
        double deltat = 1e-3;

        Eigen::MatrixXd dC(1, 3);
        P1P2(t + deltat, dC, 0);
        N.row(i) = 1./deltat * (dC - C.row(i));
        N.row(i).normalize();
    }
}

void CatmullRomChain(const Eigen::MatrixXd & p,
                     const int nsample,
                     Eigen::MatrixXd & C,
                     Eigen::MatrixXd & N) {
    if (p.rows() < 4) {
        std::cout << "not enough control points" << std::endl;
        exit(1);
    }

    // Need to determine the total number of rows and nPoints
    // for each interval first, then fill in a block of matrix
    // with the correpsonding entries

    int segnum = p.rows() -3;
    double totallength = 0;
    double curvlen[segnum];
    for (int i = 0; i < segnum; i++) {
        curvlen[i] = (p.row(i+2) - p.row(i+1)).norm();
        totallength += curvlen[i];
    }

    double avgcurvlen = totallength / nsample;
    int nPoints[segnum];
    int cumPoints[segnum];
    int totalPts = 0;
    for (int i = 0; i < segnum; i++) {
        nPoints[i] = ceil(curvlen[i] / avgcurvlen);
        cumPoints[i] = totalPts;
        totalPts += nPoints[i];
    }

    C.resize(totalPts, 3);
    N.resize(totalPts, 3);
    for (int i = 0; i < segnum; i++) {
        Eigen::MatrixXd Ci, Ni;
        int i0 = i;
        int i1 = i+1;
        int i2 = i+2;
        int i3 = i+3;

        Eigen::MatrixXd P0 = p.row(i0);
        Eigen::MatrixXd P1 = p.row(i1);
        Eigen::MatrixXd P2 = p.row(i2);
        Eigen::MatrixXd P3 = p.row(i3);
        CatmullRomSpline(P0, P1, P2, P3, Ci, Ni, nPoints[i]);
        C.block(cumPoints[i],  0, nPoints[i], 3) = Ci;
        N.block(cumPoints[i],  0, nPoints[i], 3) = Ni;
    }
}