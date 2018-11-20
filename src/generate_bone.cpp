#include "generate_bone.h"
#include <igl/read_triangle_mesh.h>
#include <math.h>
#include <Eigen/Geometry>

void generate_bones(const Eigen::MatrixXd & points,
                    std::vector<Eigen::MatrixXd> & VV,
                    std::vector<Eigen::MatrixXi> & FF) {
  
  // Create bones p(i) and p(i-1)
  for (int i = 1; i < points.rows(); i++) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //generate_flinstone_bone(points.row(i - 1), points.row(i), V, F);
    VV.push_back(V);
    FF.push_back(F);
  }
}

void generate_bone(const Eigen::Vector3d &p0,
                   const Eigen::Vector3d &p1,
                   Eigen::MatrixXd &V,
                   Eigen::MatrixXi &F) {
  // Magic number
  int part = 10;

  Eigen::Vector3d vec = (p1 - p0) / part;
  Eigen::Vector3d u, v;
  u = Eigen::Vector3d(- vec(1) / vec(0), 1, 0).normalized();
  v = u.cross(vec).normalized();

  V.resize(4 * (part + 1), 3);
  Eigen::Vector3d points[] = {p0, p0 + u, p0 + u + v, p0 + v};

  // Generating vertices
  for (int i = 0; i < part + 1; i++) {
    for (int j = 0; j < 4; j++) {
      V.row(4 * i + j) = points[j] + i * vec;
    }
  }

  F.resize(part * 8 + 4, 3);
  // Generating faces
  for (int i = 0; i < part; i++) {
    for (int j = 0; j < 4; j++) {
      int x = j;
      int y = (j + 1) % 4;
      F.row(8 * i + 2 * j + 0) = Eigen::Vector3i(4 * (i + 1) + x, 4 * i + y, 4 * i + x);
      F.row(8 * i + 2 * j + 1) = Eigen::Vector3i(4 * (i + 1) + x, 4 * (i + 1) + y, 4 * i + y);
    }
  }

  F.row(part * 8 + 0) = Eigen::Vector3i(0, 1, 2);
  F.row(part * 8 + 1) = Eigen::Vector3i(2, 3, 0);
  F.row(part * 8 + 2) = Eigen::Vector3i(4 * part + 2, 4 * part + 1, 4 * part + 0);
  F.row(part * 8 + 3) = Eigen::Vector3i(4 * part + 0, 4 * part + 3, 4 * part + 2);
}

void transform_flinstone_bone(const Eigen::Vector3d &p0,
                   const Eigen::Vector3d &p1,
                   Eigen::MatrixXd &V) {
  // It's a small mesh, so I'm just going to scale it up a bit to start
  Eigen::MatrixXd scale = Eigen::MatrixXd::Identity(3, 3) * 2.0;
  V = V * scale;
  // Find the lengths of the bone in x-y
  double bone_length_x = V.col(0).maxCoeff() - V.col(0).minCoeff();
  double bone_length_y = V.col(1).maxCoeff() - V.col(1).minCoeff();
  // Center bone better first;
  for (int i = 0; i < V.rows(); i++) { //translate
    V(i, 1) -= bone_length_y/2.0;
  }
  // Find along length (x)
  Eigen::Vector3d vec = p1 - p0;
  double stretch_x_coeff = vec.norm() / bone_length_x;
  if (stretch_x_coeff < 0.5) {
    stretch_x_coeff = 0.5;
  }
  // Find rotation
  double angle = -std::atan2(vec(1), vec(0)); // Angle btwn vector and the x-axis
  // Find translation
  Eigen::Vector3d midpoint = p0 + vec / 2.0;
  // Apply transforms
  V.col(0) = V.col(0) * stretch_x_coeff; //stretch

  Eigen::MatrixXd rotate(3, 3);
  rotate << std::cos(angle), -std::sin(angle), 0,
              std::sin(angle), std::cos(angle), 0,
              0, 0, 1;
  V = V * rotate; //rotate
  
  for (int i = 0; i < V.rows(); i++) { //translate
    V.row(i) += midpoint;
  }
}

void start_flinstone_bone(const Eigen::Vector3d &p,
                          Eigen::MatrixXd &V) {
  // It's a small mesh, so I'm just going to scale it up a bit to start
  Eigen::MatrixXd scale = Eigen::MatrixXd::Identity(3, 3) * 2.0;
  V = V * scale;
  // Find the length of the bone
  double bone_length_x = V.col(0).maxCoeff() - V.col(0).minCoeff();
  double bone_length_y = V.col(1).maxCoeff() - V.col(1).minCoeff();
  // Find translation
  Eigen::Vector3d t = p + Eigen::Vector3d(0, -bone_length_y/2.0, 0);
  for (int i = 0; i < V.rows(); i++) { //translate
    V.row(i) += t;
  }
}