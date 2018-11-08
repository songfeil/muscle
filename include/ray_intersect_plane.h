#ifndef INTRODUCTION_RAY_INTERSECT_PLANE_H
#define INTRODUCTION_RAY_INTERSECT_PLANE_H

#include <Eigen/Core>

//Find 
//Input:
//  p0      point on plane
//  n       normal of plane
//  source  of ray
//  dir     of ray
//Output
//  p   intersection point
void ray_intersect_plane(const Eigen::Vector3d &p0, const Eigen::Vector3d &n,
                         const Eigen::Vector3d &source, const Eigen::Vector3d &dir,
                         Eigen::Vector3d &p);

#endif //INTRODUCTION_GENERATE_BONE_H
