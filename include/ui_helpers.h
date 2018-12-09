

#ifndef INTRODUCTION_UI_HELPERS_H
#define INTRODUCTION_UI_HELPERS_H

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <vector>

void intersection_with_xy_plane( const igl::opengl::glfw::Viewer & viewer,
                                 const Eigen::Vector3d & last_mouse,
                                 Eigen::Vector3d & intersection );

void add_face_to_patch( const int fid,
                        const Eigen::Vector3i & face,
                        std::vector<Eigen::MatrixXi> & patches,
                        std::vector<std::vector<int>> & patch_fids );

void verts_within_x_range(const double xmin,
                          const double xmax,
                          const Eigen::MatrixXd & V,
                          std::vector<int> & vids);

#endif //INTRODUCTION_UI_HELPERS_H