//
// Created by jnyao on 11/7/18.
//

#include "ray_intersect_plane.h"

void ray_intersect_plane(const Eigen::Vector3d &p0,
                         const Eigen::Vector3d &n,
                         const Eigen::Vector3d &source,
                         const Eigen::Vector3d &dir,
                         Eigen::Vector3d &p) {
    double t = ((p0 - source).dot(n)) / dir.dot(n);
    p = source + t * dir;
}