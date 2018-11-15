// Fit n points to a flat-topped gaussian. Pass in peak (max val), plateau (corresponds to width of flat-top), and stddev (width of bell)

#ifndef INTRODUCTION_VOLUME_ALONG_CURVE_H
#define INTRODUCTION_VOLUME_ALONG_CURVE_H

#include "Eigen/Core"

void volume_along_curve(const Eigen::VectorXd curve,
                        Eigen::MatrixXd volume);

#endif //INTRODUCTION_VOLUME_ALONG_CURVE_H
