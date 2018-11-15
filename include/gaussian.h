// Fit n points to a flat-topped gaussian. Pass in peak (max val), plateau (corresponds to width of flat-top), and stddev (width of bell)

#ifndef INTRODUCTION_GAUSSIAN_H
#define INTRODUCTION_GAUSSIAN_H

#include "Eigen/Core"

void gaussian(const int n,
              const double peak,
              const double stddev,
              const double plateau,
              Eigen::VectorXd & G);

#endif //INTRODUCTION_GAUSSIAN_H
