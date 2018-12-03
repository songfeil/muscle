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
#include <igl/per_face_normals.h>
#include <Eigen/Dense>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/doublearea.h>
#include <igl/sort.h>
#include <set>

void generate_muscle(const Eigen::MatrixXd & points,
                     const int num_points,
                     const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     const std::set<int> & selected_faces,
                     std::vector<Eigen::MatrixXd> & VV,
                     std::vector<Eigen::MatrixXi> & FF) {
    Eigen::MatrixXd Vm;
    Eigen::MatrixXi Fm;
    Eigen::MatrixXd p  = points.block(points.rows() - num_points, 0, num_points, 3);
    p.colwise().reverse();
    // Your code here to populate V and F
    Eigen::MatrixXd Bc, Nc, pV, pN;
//    bezier(p1, p2, p3, 50, Bc, Nc);
//    bezier(p, num_points-1, 50, Bc, Nc);
    CatmullRomChain(p, 50, Bc, Nc);
//            bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, Bc, Nc);
    volume_along_curve(Bc, Nc, pV, pN);
    Eigen::MatrixXd All(pV.rows(), 6);
    All << pV, pN;
    poisson_surface_reconstruction(pV, pN, Vm, Fm);

    Eigen::VectorXi bb = Eigen::VectorXi(selected_faces.size()*3);
    Eigen::MatrixXd Bonedest(selected_faces.size()*3, 3);
    Eigen::MatrixXd musclef(selected_faces.size()*3, 3);

    for (auto it = selected_faces.begin(); it != selected_faces.end(); ++it) {
        int curr = std::distance(selected_faces.begin(), it);
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
        Bonedest.row(3*curr) = P.row(2);
        Bonedest.row(3*curr+1) = P.row(1);
        Bonedest.row(3*curr+2) = P.row(0);

        double minlensum = std::numeric_limits<double>::infinity();
        // Determine the optimal correspondence of vertices
        for (int i = 0; i < 3; i++) {
            double lensum = (P.row(2) - Vm.row(Fm(fi, i))).norm()
                            +(P.row(1) - Vm.row(Fm(fi, (i+1)%3 ))).norm()
                            +(P.row(0) - Vm.row(Fm(fi, (i+2)%3 ))).norm();
            if (lensum < minlensum) {
                minlensum = lensum;
                bb(3*curr) = Fm(fi, i);
                bb(3*curr+1) = Fm(fi, (i+1)%3);
                bb(3*curr+2) = Fm(fi, (i+2)%3);
                musclef.row(3*curr) = Vm.row(bb(3*curr));
                musclef.row(3*curr+1) = Vm.row(bb(3*curr+1));
                musclef.row(3*curr+2) = Vm.row(bb(3*curr+2));

            }
        }

//        bb(3*curr) = Fm(fi, 0);
//        bb(3*curr+1) = Fm(fi, 1);
//        bb(3*curr+2) = Fm(fi, 2);
//        musclef.row(3*curr) = Vm.row(bb(3*curr));
//        musclef.row(3*curr+1) = Vm.row(bb(3*curr+1));
//        musclef.row(3*curr+2) = Vm.row(bb(3*curr+2));
    }
    igl::ARAPData data;
    igl::arap_precomputation(Vm, Fm, 3, bb, data);
    //int ninterp = 10;
    int ninterp = 10;
    for (int ni = 1; ni <= ninterp; ni++) {
        double t = double (ni) / double (ninterp);
        //std::cout<< (1-t)*Fcc + t*Bcc <<std::endl;
        Eigen::MatrixXd currdest = (1-t)*musclef + t*Bonedest;
        igl::arap_solve(currdest, data, Vm);
    }
//    igl::arap_solve(Bcc, data, Vm);
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
    VV.push_back(Vm);
    FF.push_back(Fm);
    std::cout << "generate muscle" << std::endl;
}

void generate_muscle_multiface(const Eigen::MatrixXd & points,
                               const int num_points,
                               const Eigen::MatrixXd & V,
                               const Eigen::MatrixXi & F,
                               const std::set<int> & selected_faces1,
                               const std::set<int> & selected_faces2,
                               std::vector<Eigen::MatrixXd> & VV,
                               std::vector<Eigen::MatrixXi> & FF) {
    Eigen::MatrixXd Vm;
    Eigen::MatrixXi Fm;
    Eigen::MatrixXd p  = points.block(points.rows() - num_points, 0, num_points, 3);
    p.colwise().reverse();
    // Your code here to populate V and F
    Eigen::MatrixXd Bc, Nc, pV, pN;
//    bezier(p1, p2, p3, 50, Bc, Nc);
//    bezier(p, num_points-1, 50, Bc, Nc);
    CatmullRomChain(p, 50, Bc, Nc);
//            bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, Bc, Nc);
    volume_along_curve(Bc, Nc, pV, pN);
    Eigen::MatrixXd All(pV.rows(), 6);
    All << pV, pN;
    poisson_surface_reconstruction(pV, pN, Vm, Fm);

    // Calculate area and center of attachment site 1 on bone
    double area1 = 0; double num1 = 0;
    Eigen::RowVector3d center1 = Eigen::RowVector3d::Zero();
    for (auto it = selected_faces1.begin(); it != selected_faces1.end(); ++it) {
        // currently only one face is selected in face 1
        ++it;
        if (it == selected_faces1.end()) {
            break;
        }
        --it;
        // delete the above portion when selected_faces1 and selected_faces2 are correctly implemented

        num1 += 1;
        int curr = std::distance(selected_faces1.begin(), it);
        Eigen::RowVectorXi triangle = F.row(*it);
        Eigen::RowVector3d triangleCenter = (V.row(triangle(0)) + V.row(triangle(1)) + V.row(triangle(2))) / 3.0;
        center1 += triangleCenter;

        double a = (V.row(triangle(0)) - V.row(triangle(1))).norm();
        double b = (V.row(triangle(1)) - V.row(triangle(2))).norm();
        double c = (V.row(triangle(2)) - V.row(triangle(0))).norm();
        double s = (a+b+c) / 2.;
        area1 += pow(s*(s-a)*(s-b)*(s-c), 0.5);

    }
    center1 = 1./num1 * center1;

    std::cout << "num1 " << num1 <<std::endl;
    std::cout << "center1 " << center1 <<std::endl;

    // Calculate average surface area of muscle mesh
    Eigen::MatrixXd dblA;
    igl::doublearea(Vm, Fm, dblA);
    double muscle_totalarea = 0;
    for (int i = 0; i < dblA.rows(); i++) {
        muscle_totalarea += dblA(i);
    }
    double avg_area = muscle_totalarea / dblA.rows();

    // Determine the total number of faces on the muscle to be deformed towards bone
    int sourcenum1 = ceil(area1 / avg_area);
//    sourcenum1 = 100;
    std::cout << "sourcenum1 " << sourcenum1 << std::endl;
    Eigen::MatrixXd lst_dist(Fm.rows(), 1);
    triangle_hunt_lst(center1, Vm, Fm, lst_dist);
    Eigen::MatrixXd sY, sIX;
    igl::sort(lst_dist, 1, true, sY, sIX);

    Eigen::MatrixXd closest1 = sIX.block(0, 0, sourcenum1, 1);
    Eigen::MatrixXi Fpatch1(closest1.rows(), 3);
    int Vidx[3*closest1.rows()];
    int numdistinct = 0;

    // Get the total number of distinct vertices
    std::set<int> s1;
    for (int i = 0; i < closest1.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int fj = Fm(closest1(i), j);
            auto result = s1.insert(fj);
            if ( result.second ) {
                Vidx[numdistinct] = fj;
                Fpatch1(i, j) = numdistinct;
                numdistinct++;
            } else {
                for (int k = 0; k < numdistinct; k++) {
                    if (Vidx[k] == fj) {
                        Fpatch1(i, j) = k;
                    }
                }
            }
        }
    }
    Eigen::MatrixXd Vpatch1(numdistinct, 3);
    for (int i = 0; i < numdistinct; i++) {
        Vpatch1.row(i) = Vm.row(Vidx[i]);
    }

    for (int i = 0; i < closest1.rows(); i++) {
        Fpatch1.row(i) = Fm.row(closest1(i));
    }

    VV.push_back(Vm);
    FF.push_back(Fpatch1);


//    Eigen::VectorXi bnd;
//    igl::boundary_loop(F,bnd);
//
//    Eigen::MatrixXd bnd_uv;
//    igl::map_vertices_to_circle(V,bnd,bnd_uv);
//
//    Eigen::MatrixXd V_uv;
//    igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

}