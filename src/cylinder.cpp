//
// Created by 宋飞龙 on 2018-11-11.
//

#include "cylinder.h"
#include <iostream>
#include <igl/cylinder.h>

void cylinder(
    int div,
    int part,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F
) {
  div = 6;
  part = 4;

  double h = 5.0;

  // Cylinder
  Eigen::MatrixXd LV;
  Eigen::MatrixXi LF;
  igl::cylinder(30, 20, LV, LF);


//  CappedCylinderMesh cm(1.0, h, 32, part, 0.0, gml::radians(360.0));

//  LF.resize(count(cm.triangles()), 3);
//  int i = 0;
//  for (const auto& f : cm.triangles()) {
//    auto face = f.vertices.data();
//    LF.row(i) = Eigen::RowVector3i(face[0], face[1], face[2]);
//    i++;
//  }
//
//  LV.resize(count(cm.vertices()), 3);
//  i = 0;
//  for (const auto& v : cm.vertices()) {
//    auto pos = v.position.data();
//    LV.row(i) = Eigen::RowVector3d(pos[0], pos[1], pos[2]);
//    i++;
//  }
//
  Eigen::Matrix3d rot = Eigen::Matrix3d(3, 3);
  rot << 1, 0, 0,
      0, cos(0.5 * M_PI), sin(0.5 * M_PI),
      0, -sin(0.5 * M_PI), cos(0.5 * M_PI);

  LV = LV * rot.transpose();
  for (int i = 0; i < LV.rows(); i++) {
    // Scale, Rotate and translate.
    LV.row(i) = LV.row(i) + Eigen::RowVector3d(0, h / 2, 0);
  }

  V = LV;
  F = LF;

}