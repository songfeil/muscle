//
// Created by 宋飞龙 on 2018-11-11.
//

#include "pick_constrain_point.h"
#include <limits>
#include <vector>
#include <iostream>
using namespace std;

void pick_constrain_point(
    const Eigen::MatrixXd & V,
    std::vector<int> & b
) {

    double maxz = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < V.rows(); i++) {
      if (V(i, 2) > maxz) maxz = V(i, 2);
    }

  for (int i = 0; i < V.rows(); i++) {
    if (abs(V(i, 2) - maxz) < 1e-7) {
      b.emplace_back(i);
    }
  }
}