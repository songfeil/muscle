#ifndef INTRODUCTION_INFLATE_MUSCLE_H
#define INTRODUCTION_INFLATE_MUSCLE_H

#include <Eigen/Core>
#include <set>
#include <vector>

void xflate_muscle(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const int submesh_start,
                    const int submesh_end,
                    const int mode, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew);

void xflate_verts(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<int> vids,
                    const int mode, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew);

void xflate_verts_in_xrange(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const double xmid,
                    const double xrange,
                    const int mode, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew);

#endif //INTRODUCTION_INFLATE_MUSCLE_H