//
// Created by 宋飞龙 on 2018-11-11.
//

#include <iostream>
#include <igl/arap.h>
#include "generate_muscle.h"
#include "bezier.h"
#include "volume_along_curve.h"
#include "poisson_surface_reconstruction.h"
#include "triangle_hunt.h"
#include "deform.h"

void generate_muscle(const Eigen::MatrixXd & points,
                     const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     const std::set<int> & selected_faces,
                     std::vector<Eigen::MatrixXd> & VV,
                     std::vector<Eigen::MatrixXi> & FF) {

  Eigen::MatrixXd Vm;
  Eigen::MatrixXi Fm;
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
  poisson_surface_reconstruction(pV, pN, Vm, Fm);

  Eigen::VectorXi bb = Eigen::VectorXi(6);
  Eigen::MatrixXd Bcc(6, 3);

  for (auto it = selected_faces.begin(); it != selected_faces.end(); ++it) {

    Eigen::RowVectorXi triangle = F.row(*it);
    std::cout<<triangle<<std::endl;
    Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
    P.row(0) = V.row(triangle(0));
    P.row(1) = V.row(triangle(1));
    P.row(2) = V.row(triangle(2));
    std::cout<<"selected triangle"<<std::endl;
    std::cout<<P<<std::endl;
    int fi = triangle_hunts(P, Vm, Fm);
    Eigen::RowVector3i f = Fm.row(fi);
    int curr = std::distance(selected_faces.begin(), it);
    Bcc.row(3*curr) = P.row(2);
    Bcc.row(3*curr+1) = P.row(1);
    Bcc.row(3*curr+2) = P.row(0);
    bb(3*curr) = Fm(fi, 2);
    bb(3*curr+1) = Fm(fi, 1);
    bb(3*curr+2) = Fm(fi, 0);
  }

  igl::ARAPData data;
  igl::arap_precomputation(Vm, Fm, 3, bb, data);
  igl::arap_solve(Bcc, data, Vm);
  //Smooth the surface
          //    Eigen::MatrixXd Vcpy(V);
          //    Eigen::SparseMatrix<double> L, M;
          //    igl::cotmatrix(V, F, L);
          //    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
          //    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(M);
          //    Eigen::SparseMatrix<double> MinvL = solver.solve(L);
          //    Eigen::SparseMatrix<double> QL = L.transpose()*MinvL;
          //    const double al = 8e-2;
          //    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> lapSolver(al*QL + (1.-al)*M);
          //    V = lapSolver.solve(al*M*Vcpy);
          //    std::cout << V << std::endl;
          //  deform(p1, p2, p3, V, F);
  VV.push_back(Vm);
  FF.push_back(Fm);
  std::cout << "generate muscle" << std::endl;
}
