#include "volume_along_curve.h"
#include "gaussian.h"
#include <Eigen/Geometry>
#include <iostream>
#include <functional>

void generate_circle(const int n, const double radius, Eigen::MatrixXd & V) {
    V = Eigen::MatrixXd::Zero(n, 2);
    for (int i = 0; i < n; i++) {
        double theta = ((double)i / n) * 2 * M_PI;
        V.row(i) = Eigen::Vector2d(radius * cos(theta), radius * sin(theta));
    }
}

void generate_disk(const int n, const double a, const double b, Eigen::MatrixXd & V, Eigen::MatrixXd & N) {
    V = Eigen::MatrixXd::Zero(n, 2);
    N = Eigen::MatrixXd::Zero(n, 2);
    for (int i = 0; i < n; i++) {
        double theta = ((double)i / n) * 2 * M_PI;
        double r = (a * b) / sqrt(a * a * sin(theta) * sin(theta) + b * b * cos(theta) * cos(theta));
        Eigen::Vector2d pos = Eigen::Vector2d(r * cos(theta), r * sin(theta));
        V.row(i) = pos;

        // Calculate normal numerically
        double epsilon = 1e-7;
        theta += epsilon;
        double re = (a * b) / sqrt(a * a * sin(theta) * sin(theta) + b * b * cos(theta) * cos(theta));

        Eigen::Vector2d tangent = Eigen::Vector2d(r * cos(theta), r * sin(theta)) - pos;
        Eigen::Vector3d tan3d = Eigen::Vector3d(tangent(0), tangent(1), 0);
        Eigen::Vector3d znorm = Eigen::Vector3d(0, 0, 1);
        Eigen::Vector3d norm3d = znorm.cross(tan3d);
        N.row(i) = - Eigen::Vector2d(norm3d(0), norm3d(1));
    }
}

void tangent(const Eigen::Vector3d & N, Eigen::Vector3d & T, Eigen::Vector3d & B) {
    Eigen::Vector3d nn = N.normalized();
    Eigen::Vector3d x = Eigen::Vector3d(1, 0, 0);
    Eigen::Vector3d z = Eigen::Vector3d(0, 0, 1);

    if (abs(nn(0)) < 0.1 && abs(nn(1)) < 0.1) {
        T = (nn.cross(x)).normalized();
    } else {
        T = (nn.cross(z)).normalized();
    }
    B = (nn.cross(T)).normalized();
}

void volume_along_curve(const Eigen::MatrixXd & curve,
                        const Eigen::MatrixXd & normal,
                        Eigen::MatrixXd & volume,
                        Eigen::MatrixXd & pointNormal){

    int n = curve.rows();
    Eigen::VectorXd radii;
    // Fit radii of volume along gaussian distribution
    // Let's say that at it's widest, the radius should be a quarter as large as the curve is long
    // Find approximate length of the curve:
    double approx_length = 0;
    for (int i = 1; i < curve.rows(); i++) {
        approx_length += (curve.row(i) - curve.row(i - 1)).norm();
    }
    gaussian(n, approx_length / 4.0, approx_length, 2.0, radii);
    // Next we need to sample points in rings of radius r along the besier
    // Option 1) create a function for a circle on a plane perpendicular to the tangent and centered at each bezier-point with its corresopnding radius, sample along that
    // Option 2) find a vector of length r perpendicular to the tangent at each bezier-point, rotate it 360 degrees.
    // Either way we need to find normals to points on the curve
    // For now, our curves only exist on the X-Y plane (z = 0 for all points) just based on the user interface
    // So we can use that to our advantage, and generalize it later...
    int circleSampleCount = 15;

    auto radiusFunc = [](double max, double t){
        return max * (-3 * (t - 0.5) * (t - 0.5) + 1);
    };

    volume = Eigen::MatrixXd::Zero(n * circleSampleCount, 3);
    pointNormal = Eigen::MatrixXd::Zero(n * circleSampleCount, 3);
    for (int i = 0; i < n; i++) {
        double radius = radiusFunc(0.5, (double)i/n);
        Eigen::MatrixXd circle2D, circle2dNorm;
//      generate_circle(circleSampleCount, radius, circle2D);
        generate_disk(circleSampleCount, radius, radius * 0.5, circle2D, circle2dNorm);

        Eigen::Vector3d N, T, B, pos;
        N = normal.row(i);
        N.normalize();
        pos = normal.row(i);

        tangent(N, T, B);

        Eigen::MatrixXd transform = Eigen::MatrixXd(3, 2);
        transform.col(0) = T;
        transform.col(1) = B;
        std::cout << "Normal " << N << std::endl;
        std::cout << transform << std::endl << std::endl;

        Eigen::MatrixXd scaleCircle2D = circle2D;
        Eigen::MatrixXd tCircle = (transform * scaleCircle2D.transpose()).transpose();
        Eigen::MatrixXd tNorm = (transform * circle2dNorm.transpose()).transpose();
//      Eigen::MatrixXd tCircle = Eigen::MatrixXd(circle2D.rows(), 3);
//      for (int j = 0; j < tCircle.rows(); j++) {
//        Eigen::Vector3d pt = circle2D.row(j);
//        tCircle.row(j) = transform * pt;
//      }
        for (int j = 0; j < tCircle.rows(); j++) {
            tCircle.row(j) += curve.row(i);
            volume.row(i * circleSampleCount + j) = tCircle.row(j);
//        pointNormal.row(i * circleSampleCount + j) = volume.row(i * circleSampleCount + j) - curve.row(i);
            pointNormal.row(i * circleSampleCount + j) = tNorm.row(j);
        }
    }
}

void ellipse_along_curve(const Eigen::MatrixXd & curve,
                         const Eigen::MatrixXd & normal,
                         const Eigen::VectorXd & long_axis,
                         const Eigen::MatrixXd & long_dir,
                         const Eigen::VectorXd & short_axis,
                         const Eigen::MatrixXd & short_dir,
                         Eigen::MatrixXd & volume,
                         Eigen::MatrixXd & pointNormal){

    int n = curve.rows();
    Eigen::VectorXd radii;
    // Fit radii of volume along gaussian distribution
    // Let's say that at it's widest, the radius should be a quarter as large as the curve is long
    // Find approximate length of the curve:
    double approx_length = 0;
    for (int i = 1; i < curve.rows(); i++) {
        approx_length += (curve.row(i) - curve.row(i - 1)).norm();
    }
    gaussian(n, approx_length / 4.0, approx_length, 2.0, radii);
    // Next we need to sample points in rings of radius r along the besier
    // Option 1) create a function for a circle on a plane perpendicular to the tangent and centered at each bezier-point with its corresopnding radius, sample along that
    // Option 2) find a vector of length r perpendicular to the tangent at each bezier-point, rotate it 360 degrees.
    // Either way we need to find normals to points on the curve
    // For now, our curves only exist on the X-Y plane (z = 0 for all points) just based on the user interface
    // So we can use that to our advantage, and generalize it later...
    int circleSampleCount = 15;

    // Define a function that reach the minimum at the midpoint
    auto radiusFunc = [](double max, double t){
        return max * (-3 * (t - 0.5) * (t - 0.5) + 1);
    };

    volume = Eigen::MatrixXd::Zero(n * circleSampleCount, 3);
    pointNormal = Eigen::MatrixXd::Zero(n * circleSampleCount, 3);
    for (int i = 0; i < n; i++) {
        double radius = radiusFunc(0.5, (double)i/n);
        Eigen::MatrixXd circle2D, circle2dNorm;
//      generate_circle(circleSampleCount, radius, circle2D);
        generate_disk(circleSampleCount, long_axis(i), short_axis(i), circle2D, circle2dNorm);

        Eigen::Vector3d N, T, B, pos;
        N = normal.row(i);
        N.normalize();
        pos = normal.row(i);

        tangent(N, T, B);

        Eigen::MatrixXd transform = Eigen::MatrixXd(3, 2);
        transform.col(0) =  long_dir.row(i).transpose(); // long axis
        transform.col(1) = short_dir.row(i).transpose(); // short axis
        std::cout << "Normal " << N << std::endl;
        std::cout << transform << std::endl << std::endl;

        Eigen::MatrixXd scaleCircle2D = circle2D;
        Eigen::MatrixXd tCircle = (transform * scaleCircle2D.transpose()).transpose();

        // Normal of the surface may not be the normal of the ellipse transformed
        Eigen::MatrixXd tNorm = (transform * circle2dNorm.transpose()).transpose();
//      Eigen::MatrixXd tCircle = Eigen::MatrixXd(circle2D.rows(), 3);
//      for (int j = 0; j < tCircle.rows(); j++) {
//        Eigen::Vector3d pt = circle2D.row(j);
//        tCircle.row(j) = transform * pt;
//      }
        for (int j = 0; j < tCircle.rows(); j++) {
            tCircle.row(j) += curve.row(i);
            volume.row(i * circleSampleCount + j) = tCircle.row(j);
//        pointNormal.row(i * circleSampleCount + j) = volume.row(i * circleSampleCount + j) - curve.row(i);
            pointNormal.row(i * circleSampleCount + j) = tNorm.row(j);
        }
    }
}