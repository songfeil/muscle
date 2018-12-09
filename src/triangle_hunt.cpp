//
// Created by 宋飞龙 on 2018-11-15.
//

#include "triangle_hunt.h"
#include <limits>
#include <iostream>
#include <math.h>

using namespace std;

void triangle_hunt(Eigen::Matrix3d & P, Eigen::MatrixXd & V, Eigen::MatrixXi & F) {
    Eigen::RowVector3d triangleCenter = (P.row(0) + P.row(1) + P.row(2)) / 3.0;
    double dist = std::numeric_limits<double>::infinity();
    int fi = -1;

    for (int i = 0; i < F.rows(); i++) {
        Eigen::RowVector3i f = F.row(i);
        Eigen::RowVector3d p0 = V.row(f(0));
        Eigen::RowVector3d p1 = V.row(f(1));
        Eigen::RowVector3d p2 = V.row(f(2));

        Eigen::RowVector3d pc = (p0 + p1 + p2) / 3.0;

        double cdist = abs((pc - triangleCenter).norm());

        if (cdist < dist) {
            dist = cdist;
            fi = i;
        }
    }

    // Change the location
    Eigen::RowVector3i f = F.row(fi);
    cout<<"triangle hunt V"<<endl;
    cout<<V.row(f(0))<<endl;
    V.row(f(0)) = P.row(0);
    V.row(f(1)) = P.row(1);
    V.row(f(2)) = P.row(2);
    cout<<"triangle hunt P"<<endl;
    cout<<P.row(0)<<endl;
}

int triangle_hunts(Eigen::Matrix3d & P, Eigen::MatrixXd V, Eigen::MatrixXi F) {
    Eigen::RowVector3d triangleCenter = (P.row(0) + P.row(1) + P.row(2)) / 3.0;
    std::cout << "triangle center" << std::endl;
    std::cout << triangleCenter << std::endl;
    double dist = std::numeric_limits<double>::infinity();
    int fi = -1;
    Eigen::RowVector3d sp0 = P.row(2);
    Eigen::RowVector3d sp1 = P.row(1);
    Eigen::RowVector3d sp2 = P.row(0);

    for (int i = 0; i < F.rows(); i++) {
        Eigen::RowVector3i f = F.row(i);

//        Eigen::RowVector3d pc = (p0 + p1 + p2) / 3.0;

        double minlensum = std::numeric_limits<double>::infinity();
        for (int k = 0; k < 3; k++) {
            Eigen::RowVector3d p0 = V.row(f(k));
            Eigen::RowVector3d p1 = V.row(f((k+1)%3));
            Eigen::RowVector3d p2 = V.row(f((k+2)%3));
            double lensum = (sp0 - p0 ).norm()
                            +(sp1 - p1 ).norm()
                            +(sp2 - p2 ).norm();
            if (lensum < minlensum) {
                minlensum = lensum;
            }
        }

//        double cdist = abs((pc - triangleCenter).norm());
        double cdist = minlensum;

        if (cdist < dist) {
            dist = cdist;
            fi = i;
        }
    }

    // Change the location
    Eigen::RowVector3i f = F.row(fi);
    return fi;
}


void triangle_hunt_lst(Eigen::RowVector3d & P, const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::MatrixXd & lst_dist) {
    // Record the distance from the center of each triangle on muscle to center point P of the selected faces on bone
    double dist = std::numeric_limits<double>::infinity();
    int fi = -1;
    for (int i = 0; i < F.rows(); i++) {
        Eigen::RowVector3i f = F.row(i);
        Eigen::RowVector3d p0 = V.row(f(0));
        Eigen::RowVector3d p1 = V.row(f(1));
        Eigen::RowVector3d p2 = V.row(f(2));

        Eigen::RowVector3d pc = (p0 + p1 + p2) / 3.0;

        double cdist = abs((pc - P).norm());
        lst_dist(i, 0) = cdist;
        if (cdist < dist) {
            dist = cdist;
            fi = i;
        }
    }
}