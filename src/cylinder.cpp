//
// Created by 宋飞龙 on 2018-11-11.
//

#include "cylinder.h"
#include <generator/generator.hpp>
#include <generator/CylinderMesh.hpp>
#include <iostream>

using namespace generator;

void cylinder(
    int div,
    int part,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F
) {
  div = 6;
  part = 10;

  CylinderMesh cm;

  F.resize(count(cm.triangles()), 3);
  int i = 0;
  for (const auto& f : cm.triangles()) {
    auto face = f.vertices.data();
    F.row(i) = Eigen::RowVector3i(face[0], face[1], face[2]);
    i++;
  }

  V.resize(count(cm.vertices()), 3);
  i = 0;
  for (const auto& v : cm.vertices()) {
    auto pos = v.position.data();
    V.row(i) = Eigen::RowVector3d(pos[0], pos[1], pos[2]);
    i++;
  }
}