//
// Created by jnyao on 11/7/18.
//

#include "bezier.h"
#include <math.h>
#include <iostream>
#include <igl/cat.h>

int factorial(int n) {
    if (n == 0) {
        return 1;
    } else {
        return factorial(n - 1) * n;
    }
}

int n_choose_k(int n, int k) {
    return factorial(n) / (factorial(k) * factorial(n-k));
}

void CatmullRomSpline(Eigen::MatrixXd P0, Eigen::MatrixXd P1, Eigen::MatrixXd P2, Eigen::MatrixXd P3, Eigen::MatrixXd & C, Eigen::MatrixXd & N,  int nPoints=10) {
    std::cout << "Enter CatmullRomSpline" << std::endl;
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

//        std::cout << "A" << std::endl;
//        std::cout << A1 << std::endl;
//        std::cout << A2 << std::endl;
//        std::cout << A3 << std::endl;
//        std::cout << "B" << std::endl;
//        std::cout << B1 << std::endl;
//        std::cout << B2 << std::endl;
//        std::cout << "C" << std::endl;
        Cc.row(rowi) = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2;
    };

    for (int i = 0; i < nPoints; i++) {
        double t = t1 + (t2 - t1) / (nPoints) * i + (t2 - t1) / (2*nPoints);

        P1P2(t, C, i);
        double deltat = 1e-3;

        Eigen::MatrixXd dC(1, 3);
        P1P2(t + deltat, dC, 0);
        N.row(i) = 1./deltat * (dC - C.row(i));
//        std::cout << N.row(i) << std::endl;
        N.row(i).normalize();
    }
    std::cout << "C" << std::endl;
    for (int i = 0; i < nPoints; i++) {
        std::cout << C.row(i) << std::endl;
    }
    std::cout << "N" << std::endl;
    for (int i = 0; i < nPoints; i++) {
        std::cout << N.row(i) << std::endl;
    }
}

void CatmullRomChain(const Eigen::MatrixXd & p,
//            const Eigen::Vector3d & p1,
//            const Eigen::Vector3d & p2,
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
    // TODO determine nPoints adaptively
    for (int i = 0; i < segnum; i++) {
        Eigen::MatrixXd Ci, Ni;
        int i0 = i;
        int i1 = i+1;
        int i2 = i+2;
        int i3 = i+3;

//        if (i == 0) {
//            i0 = 0; i1 = 0; i2 = 1; i3 = 2;
//        } else if (i == p.rows() - 1) {
//            i0 = p.rows() - 3; i1 = p.rows() - 2; i2 = p.rows() - 1; i3 = p.rows() - 1;
//        }

        Eigen::MatrixXd P0 = p.row(i0);
        Eigen::MatrixXd P1 = p.row(i1);
        Eigen::MatrixXd P2 = p.row(i2);
        Eigen::MatrixXd P3 = p.row(i3);
        CatmullRomSpline(P0, P1, P2, P3, Ci, Ni, nPoints[i]);
        C.block(cumPoints[i],  0, nPoints[i], 3) = Ci;
        N.block(cumPoints[i],  0, nPoints[i], 3) = Ni;
//        igl::cat(1, C, Ci, C);
//        igl::cat(1, N, Ni, N);
    }
    std::cout << "C cat" << std::endl;
    for (int i = 0; i < C.rows(); i++) {
        std::cout << C.row(i) << std::endl;
    }
    std::cout << "N cat" << std::endl;
    for (int i = 0; i < N.rows(); i++) {
        std::cout << N.row(i) << std::endl;
    }
}

void bezier(const Eigen::MatrixXd & p,
//            const Eigen::Vector3d & p1,
//            const Eigen::Vector3d & p2,
            const int ndegree,
            const int nsample,
            Eigen::MatrixXd & B,
            Eigen::MatrixXd & N) {
    Eigen::VectorXd bezier_coeff(ndegree + 1);
    for (int i = 0; i <= ndegree; i++) {
        std::cout << n_choose_k(ndegree, i) <<" "<<std::endl;
        bezier_coeff(i) = n_choose_k(ndegree, i);
    }
    auto bfunc = [&](double t){
        Eigen::RowVector3d ret_p = Eigen::Vector3d::Zero();
        for (int i = 0; i <= ndegree; i++) {
            ret_p = ret_p + pow((1-t), ndegree - i) * pow(t, i)*bezier_coeff(i) * p.row(i);
        }
//        return (1 - t) * (1 - t) * p0 + 2 * t * (1 - t) * p1 + t * t * p2;
        return ret_p;
    };
    auto nfunc = [&](double t){
        // to avoid numerical issue, t should not be 0 or 1
        Eigen::RowVector3d ret_p = Eigen::Vector3d::Zero();
        for (int i = 0; i <= ndegree; i++) {
            ret_p = ret_p + ((i-ndegree)*pow((1-t), ndegree - i - 1) * pow(t,   i) +
                                       i*pow((1-t), ndegree - i)     * pow(t, i-1)) *
                                      bezier_coeff(i) * p.row(i);
        }
//        return (2 * t - 2) * p0 + (2 - 4 * t) * p1 + 2 * t * p2;
        return ret_p;
    };

    B.resize(nsample + 1, 3);
    N.resize(nsample + 1, 3);
    for (int i = 0; i < nsample + 1; i++ ) {
        B.row(i) = bfunc(0.1 + 0.8 * ((double) i / nsample));
        N.row(i) = nfunc(0.1 + 0.8 * ((double) i / nsample));
    }

}