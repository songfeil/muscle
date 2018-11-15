#include "gaussian.h"
#include <math.h>
#include <iostream>

void gaussian(const int n,
              const double peak,
              const double width,
              const double plateau,
              Eigen::VectorXd & G) {

    G.resize(n);
    double mid = (n - 1.0) / 2.0;
    std::cout << "gaussian" << std::endl;
    for (int x = 0; x < n; x++) {
        G(x) = peak * std::exp((-1.0) * std::pow((x - mid)/width, 2 * plateau));
        std::cout << G(x) << std::endl;
    }

}