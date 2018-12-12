#ifndef INTRODUCTION_MESH_EDITING_H
#define INTRODUCTION_MESH_EDITING_H

#include <Eigen/Core>
#include <set>
#include <vector>

void xflate_mesh(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const int submesh_start,
                    const int submesh_end,
                    const int dir, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew);

void xflate_verts(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<int> vids,
                    const int dir, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew);

void xflate_verts_in_xrange(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const double xmid,
                    const double xrange,
                    const int dir, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew);

void smooth_mesh_with_fixed(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::set<int> & attached_vids,
                    Eigen::MatrixXd & Vnew);

#endif //INTRODUCTION_MESH_EDITING_H