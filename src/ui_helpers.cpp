#include "ui_helpers.h"
#include <igl/unproject_ray.h>

void prune_input_stroke(const double thresh,
                        Eigen::MatrixXd & points) {
    Eigen::MatrixXd pruned(points.rows(), 3);
    pruned.row(0) = points.row(0);
    Eigen::Vector3d point1 = points.row(0);
    Eigen::Vector3d point2;
    int i = 1;
    int n = 1;
    while (i < points.rows()) {
        point2 = points.row(i);
        if ((point1 - point2).norm() > thresh) {
            pruned.row(n) = point2;
            point1 = point2;
            n++;
        }
        i++;
    }
    pruned.conservativeResize(n + 1, 3);
    points = pruned;

}

void ray_intersect_plane(const Eigen::Vector3d &p0,
                         const Eigen::Vector3d &n,
                         const Eigen::Vector3d &source,
                         const Eigen::Vector3d &dir,
                         Eigen::Vector3d &p) {
    double t = ((p0 - source).dot(n)) / dir.dot(n);
    p = source + t * dir;
}

void intersection_with_xy_plane( const igl::opengl::glfw::Viewer & viewer,
                                 const Eigen::Vector3d & last_mouse,
                                 Eigen::Vector3d & intersection ) {
  Eigen::Vector3d source, dir;
  Eigen::Vector2d p;
  p << last_mouse.head(2)(0), last_mouse.head(2)(1);
  igl::unproject_ray(p, viewer.core.view, viewer.core.proj, viewer.core.viewport, source, dir);
  ray_intersect_plane(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), source, dir, intersection);
}

void add_face_to_patch( const int fid,
                        const Eigen::Vector3i & face,
                        std::vector<Eigen::MatrixXi> & patches,
                        std::vector<std::vector<int>> & patch_fids) {

    if (patches.size() > 0) {

        // If there are any existing patches to compare with,
        // Check if this face should be added to an existing patch

        int patch_num = 0;
        for (auto & patch: patches) {
            // Iterate over patches
            // Gather their unique verts
            std::set<int> patch_verts;
            for (int i = 0; i < patch.rows(); i++) {
                for (int j = 0; j < 3; j++) {
                    patch_verts.insert(patch(i, j));
                }
            }
            // Check edge adjacency (sharing of 2 vertices) for each patch
            int num_shared_verts = 0;
            for (int vert: patch_verts) {
                for (int i = 0; i < 3; i++) {
                    if (face(i) == vert) {
                        // If we share a vert, increment
                        num_shared_verts++;
                        if (num_shared_verts > 1) {
                            // If we have 2 shared verts, we probably share an edge -> add to this patch!
                            patch.conservativeResize(patch.rows() + 1, 3);
                            patch.row(patch.rows() - 1) = face;
                            patch_fids.at(patch_num).push_back(fid);
                            return;
                        }
                    }
                }
            }
            patch_num++;
        }
    }

  // If we didn't add this face to an existing patch, we make a new patch
  Eigen::MatrixXi new_patch(1, 3);
  new_patch.row(0) = face;
  std::vector<int> new_patch_fids;
  new_patch_fids.push_back(fid);

  patches.push_back(new_patch);
  patch_fids.push_back(new_patch_fids);

}

void verts_within_x_range(const double xmin,
                          const double xmax,
                          const Eigen::MatrixXd & V,
                          std::vector<int> & vids) {
    double x;
    for (int i = 0; i < V.rows(); i++) {
        x = V(i, 0);
        
        if ((x < xmax) && (x > xmin)) {
            vids.push_back(i);
        }
    }
}