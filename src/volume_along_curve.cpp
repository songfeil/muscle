#include "volume_along_curve.h"
#include "gaussian.h"

void volume_along_curve(const Eigen::MatrixXd curve,
                        Eigen::MatrixXd volume){

    int n = curve.rows();
    VectorXd radii;
    // Fit radii of volume along gaussian distribution
    // Let's say that at it's widest, the radius should be a quarter as large as the curve is long
    // Find approximate length of the curve:
    double approx_length = 0;
    for (int i = 1; i < curve.rows(); i++) {
        approx_length += (curve(i) - curve(i - 1)).norm();
    }
    gaussian(n, approx_length / 4.0, approx_length, 2.0, radii);
    // Next we need to sample points in rings of radius r along the besier
    // Option 1) create a function for a circle on a plane perpendicular to the tangent and centered at each bezier-point with its corresopnding radius, sample along that
    // Option 2) find a vector of length r perpendicular to the tangent at each bezier-point, rotate it 360 degrees.
    // Either way we need to find normals to points on the curve
    // For now, our curves only exist on the X-Y plane (z = 0 for all points) just based on the user interface
    // So we can use that to our advantage, and generalize it later...

}