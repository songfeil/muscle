//
// Created by jnyao on 11/7/18.
//

#include "bezier.h"

void bezier(const Eigen::Vector3d & p0,
            const Eigen::Vector3d & p1,
            const Eigen::Vector3d & p2,
            const int n,
            Eigen::MatrixXd & B) {
    auto bfunc = [&](double t){
        return (1 - t) * (1 - t) * p0 + 2 * t * (1 - t) * p1 + t * t * p2;
    };

    B.resize(n + 1, 3);
    for (int i = 0; i < n + 1; i++ ) {
        B.row(i) = bfunc((double) i / n);
    }
}