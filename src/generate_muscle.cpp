//
// Created by 宋飞龙 on 2018-11-11.
//

#include <iostream>
#include "generate_muscle.h"
#include "bezier.h"
#include "volume_along_curve.h"
#include "poisson_surface_reconstruction.h"

void generate_muscle(const Eigen::MatrixXd & points,
                    std::vector<Eigen::MatrixXd> & VV,
                    std::vector<Eigen::MatrixXi> & FF) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::Vector3d p1 = points.row(points.rows() - 1);
  Eigen::Vector3d p2 = points.row(points.rows() - 2);
  Eigen::Vector3d p3 = points.row(points.rows() - 3);
  // Your code here to populate V and F
  Eigen::MatrixXd Bc, Nc, pV, pN;
  bezier(p1, p2, p3, 50, Bc, Nc);
//            bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, Bc, Nc);
  volume_along_curve(Bc, Nc, pV, pN);
  Eigen::MatrixXd All(pV.rows(), 6);
  All << pV, pN;
  poisson_surface_reconstruction(pV, pN, V, F);
  // Smooth the surface
//              Eigen::MatrixXd Vcpy(V);
//              Eigen::SparseMatrix<double> L, M;
//              igl::cotmatrix(V, F, L);
//              igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
//              Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(M);
//              Eigen::SparseMatrix<double> MinvL = solver.solve(L);
//              Eigen::SparseMatrix<double> QL = L.transpose()*MinvL;
//              const double al = 8e-2;
//              Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> lapSolver(al*QL + (1.-al)*M);
//              V = lapSolver.solve(al*M*Vcpy);
//              std::cout << V << std::endl;
//            deform(p1, p2, p3, V, F);
  VV.push_back(V);
  FF.push_back(F);
  std::cout << "generate muscle" << std::endl;
}
