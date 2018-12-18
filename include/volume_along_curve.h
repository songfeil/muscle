// Fit n points to a flat-topped gaussian. Pass in peak (max val), plateau (corresponds to width of flat-top), and stddev (width of bell)

#ifndef INTRODUCTION_VOLUME_ALONG_CURVE_H
#define INTRODUCTION_VOLUME_ALONG_CURVE_H

#include "Eigen/Core"

// Given points along a curve, genereate point clouds and normals from circles
// around the curve for the purpose of Poisson Surface Reconstruction.
void volume_along_curve(const Eigen::MatrixXd & curve,
                        const Eigen::MatrixXd & normal,
                        Eigen::MatrixXd & volume,
                        Eigen::MatrixXd & pointNormal);

// Given points and other interpolated ellipse parameters along a curve, genereate
// point clouds and normals from ellipses around the curve for the purpose of Poisson Surface Reconstruction.
void ellipse_along_curve(const Eigen::MatrixXd & curve,
                        const Eigen::MatrixXd & normal,
                        const Eigen::VectorXd & long_axis,
                        const Eigen::MatrixXd & long_dir,
                        const Eigen::VectorXd & short_axis,
                        const Eigen::MatrixXd & short_dir,
                        Eigen::MatrixXd & volume,
                        Eigen::MatrixXd & pointNormal);
/* each of the parameters are defined point-wise on the curve along which the volume is generated
 */
#endif //INTRODUCTION_VOLUME_ALONG_CURVE_H
