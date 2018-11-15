//
// Created by 宋飞龙 on 2018-11-11.
//

#ifndef INTRODUCTION_PICK_CONSTRAIN_POINT_H
#define INTRODUCTION_PICK_CONSTRAIN_POINT_H

#include <Eigen/Core>
#include <vector>

void pick_constrain_point(
        const Eigen::MatrixXd & V,
        std::vector<int> & b
    );

#endif //INTRODUCTION_PICK_CONSTRAIN_POINT_H
