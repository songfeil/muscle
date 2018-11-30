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

    auto P1P2 = [&] (double t) {
        Eigen::MatrixXd A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1;
        Eigen::MatrixXd A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2;
        Eigen::MatrixXd A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3;

        Eigen::MatrixXd B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2;
        Eigen::MatrixXd B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3;

        return (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2;

    };

    for (int i = 0; i < nPoints; i++) {
        double t = t1 + (t2 - t1) / (nPoints-1) * i;
        C.row(i) = P1P2(t);
        double deltat = 1e-4;

        Eigen::MatrixXd dC = P1P2(t + deltat);
        N.row(i) = 1./deltat * (C.row(i) - dC);
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

    C.resize(0, 0);
    N.resize(0, 0);
    // TODO determine nPoints adaptively
    for (int i = 0; i < p.rows() - 3; i++) {
        Eigen::MatrixXd Ci, Ni;
        CatmullRomSpline(p.row(i), p.row(i+1), p.row(i+2), p.row(i+3), Ci, Ni);
        igl::cat(1, C, Ci, C);
        igl::cat(1, N, Ni, N);
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