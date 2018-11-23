//
// Created by jnyao on 11/7/18.
//

#include "bezier.h"
#include <math.h>
#include <iostream>

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