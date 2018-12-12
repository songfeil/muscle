#include "mesh_editing.h"
#include <igl/per_vertex_normals.h>
#include <igl/per_vertex_attribute_smoothing.h>
#include <cmath>

 void xflate_mesh(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const int submesh_start,
                    const int submesh_end,
                    const int dir, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew) {
     // Compute desired locations (V + unit_normals!)
    Eigen::MatrixXd N, V_new;
    igl::per_vertex_normals(V, F, N);
     // Only inflate submesh
    Vnew.resizeLike(V);
    for (int i = 0; i < V.rows(); i++){
      Vnew.row(i) = V.row(i);
      if (i >= submesh_start) {
        Vnew.row(i) += N.row(i) * dir * 0.1;
      }
    }
 }

 void xflate_verts(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<int> vids,
                    const int dir, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew) {
  Eigen::MatrixXd N, V_new;
  igl::per_vertex_normals(V, F, N);
    // Only inflate submesh
  Vnew = V.replicate(1, 1);
  for (int i : vids) {
    Vnew.row(i) += N.row(i) * dir * 0.01;
  }
}

void xflate_verts_in_xrange(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const double xmid,
                    const double xrange,
                    const int dir, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew){
  Eigen::MatrixXd N, V_new;
  igl::per_vertex_normals(V, F, igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, N);
  Vnew = V.replicate(1, 1);
  for (int i = 0; i < Vnew.rows(); i++) {
    double x = Vnew(i, 0);
    double dist = std::abs(x - xmid);
    if (dist < (xrange * 2)) {
      // Multiply by gaussian for smoothness
      double coeff = std::exp((-1.0) * std::pow((x - xmid)/(xrange / 2.0), 4));
      Vnew.row(i) += N.row(i) * dir * 0.1 * coeff;
    }
  }
}

void smooth_mesh_with_fixed(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::set<int> & attached_vids,
                    Eigen::MatrixXd & Vnew){
  Vnew.resize(V.rows(), 3);

  // Gather coordinates for the fixed verts 
  Eigen::MatrixXd fixed_verts(attached_vids.size(), 3);
  for (auto it = attached_vids.begin(); it != attached_vids.end(); ++it){
    int vid = *it;
    int curr = std::distance(attached_vids.begin(), it);
    fixed_verts.row(curr) = V.row(vid);
  }

  // Smooth normally
  igl::per_vertex_attribute_smoothing(V, F, Vnew);

  // Pin fixed verts back into place
  // HACKY, but works well for our purposes
  for (auto it = attached_vids.begin(); it != attached_vids.end(); ++it){
    int vid = *it;
    int curr = std::distance(attached_vids.begin(), it);
    Vnew.row(vid) = fixed_verts.row(curr);
  }

}