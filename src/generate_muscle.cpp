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
#include <vector>
#include <move_patch.h>
#include <PCA_elipse.h>
#include <igl/per_vertex_attribute_smoothing.h>

void generate_muscle(const Eigen::MatrixXd & points,
                     const int num_points,
                     const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     const std::vector<Eigen::MatrixXi> & selected_faces,
                     Eigen::MatrixXd & Vm,
                     Eigen::MatrixXi & Fm,
                     std::set<int> & fixed_vids) {
    Eigen::MatrixXd Vm_before_smooth;
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
    poisson_surface_reconstruction(pV, pN, Vm_before_smooth, Fm);
    // Smooth before doing any of the face-matching operations
    igl::per_vertex_attribute_smoothing(Vm_before_smooth, Fm, Vm);

    Eigen::VectorXi bb = Eigen::VectorXi(selected_faces.size()*3);
    Eigen::MatrixXd Bonedest(selected_faces.size()*3, 3);
    Eigen::MatrixXd musclef(selected_faces.size()*3, 3);


    for (auto it = selected_faces.begin(); it != selected_faces.end(); ++it) {
        int curr = std::distance(selected_faces.begin(), it);
        Eigen::RowVectorXi triangle = (*it).row(0);
        std::cout<<triangle<<std::endl;
        Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
        P.row(0) = V.row(triangle(0));
        P.row(1) = V.row(triangle(1));
        P.row(2) = V.row(triangle(2));
        std::cout<<"selected triangle"<<std::endl;
        std::cout<<P<<std::endl;
        int fi = triangle_hunts(P, Vm, Fm);
        Eigen::RowVector3i f = Fm.row(fi);
        fixed_vids.insert(f(0));
        fixed_vids.insert(f(1));
        fixed_vids.insert(f(2));
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
    
}



void generate_muscle_multiface(const Eigen::MatrixXd & points,
                               const int num_points,
                               const Eigen::MatrixXd & V,
                               const Eigen::MatrixXi & F,
                               const std::vector<Eigen::MatrixXi> & selected_faces,
                               Eigen::MatrixXd & Vm,
                               Eigen::MatrixXi & Fm) {

    Eigen::MatrixXd Vm_before_smooth;

    int large_num = 1000;
    Eigen::VectorXi bb = Eigen::VectorXi(large_num);
    Eigen::MatrixXd Bonedest(large_num, 3);
    Eigen::MatrixXd musclef(large_num, 3);
    int cum_bnd = 0;
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
    poisson_surface_reconstruction(pV, pN, Vm_before_smooth, Fm);
    igl::per_vertex_attribute_smoothing(Vm_before_smooth, Fm, Vm);


    // Calculate average surface area of muscle mesh
    Eigen::MatrixXd dblA;
    igl::doublearea(Vm, Fm, dblA);
    double muscle_totalarea = 0;
    for (int i = 0; i < dblA.rows(); i++) {
        muscle_totalarea += dblA(i);
    }
    double avg_area = muscle_totalarea / dblA.rows();


    // Calculate area and center of attachment site 1 on bone
    // Also get the total number of distinct vertices in order
    // to construct a small mesh for the patch on bone.
    for (auto it = selected_faces.begin(); it != selected_faces.end(); ++it) {
        Eigen::MatrixXi curr_patch = *it;
        double area1 = 0; double num1 = 0;
        Eigen::MatrixXi Fpatch1(curr_patch.rows(), 3);
        int Vidx[3*curr_patch.rows()];
        std::set<int> s1;
        int numdistinct = 0;
        Eigen::RowVector3d center1 = Eigen::RowVector3d::Zero();
        int curr = std::distance(selected_faces.begin(), it);
//        if (curr == 1) {
//            break;
//        }
        for (int i = 0; i < curr_patch.rows(); i++) {
            num1 += 1;
            Eigen::RowVectorXi triangle = curr_patch.row(i);
            Eigen::RowVector3d triangleCenter = (V.row(triangle(0)) + V.row(triangle(1)) + V.row(triangle(2))) / 3.0;
            center1 += triangleCenter;

            double a = (V.row(triangle(0)) - V.row(triangle(1))).norm();
            double b = (V.row(triangle(1)) - V.row(triangle(2))).norm();
            double c = (V.row(triangle(2)) - V.row(triangle(0))).norm();
            double s = (a+b+c) / 2.;
            area1 += pow(s*(s-a)*(s-b)*(s-c), 0.5);

            // Construct smaller mesh for the patch on bone
            for (int j = 0; j < 3; j++) {
                int fj = triangle(j);
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
        center1 = 1./num1 * center1;
        Eigen::MatrixXd Vpatch1(numdistinct, 3);
        for (int i = 0; i < numdistinct; i++) {
            Vpatch1.row(i) = V.row(Vidx[i]);
        }

        std::cout<< "bone patch vertices\n";
        std::cout<< Vpatch1 << std::endl;

        // Determine the total number of faces on the muscle to be deformed towards bone
        int sourcenum1 = ceil(area1 / avg_area);
        sourcenum1 = 1;
        Eigen::MatrixXd lst_dist(Fm.rows(), 1);
        triangle_hunt_lst(center1, Vm, Fm, lst_dist);
        Eigen::MatrixXd sY, sIX;
        igl::sort(lst_dist, 1, true, sY, sIX);

        Eigen::MatrixXd closest1 = sIX.block(0, 0, sourcenum1, 1);
        Eigen::MatrixXi Fmpatch1(closest1.rows(), 3);
        int Vmidx[3*closest1.rows()];
        numdistinct = 0;

        // Get the total number of distinct vertices
        std::set<int> sm1;
        for (int i = 0; i < closest1.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                int fj = Fm(closest1(i), j);
                auto result = sm1.insert(fj);
                if ( result.second ) {
                    Vmidx[numdistinct] = fj;
                    Fmpatch1(i, j) = numdistinct;
                    numdistinct++;
                } else {
                    for (int k = 0; k < numdistinct; k++) {
                        if (Vmidx[k] == fj) {
                            Fmpatch1(i, j) = k;
                        }
                    }
                }
            }
        }
        Eigen::MatrixXd Vmpatch1(numdistinct, 3);
        for (int i = 0; i < numdistinct; i++) {
            Vmpatch1.row(i) = Vm.row(Vmidx[i]);
        }

        Eigen::MatrixXd musVNew = Vmpatch1;
        Eigen::RowVector3d bonCenter = Vpatch1.colwise().mean();
        Eigen::RowVector3d musCenter = Vmpatch1.colwise().mean();
        Eigen::Vector3d t = bonCenter - musCenter;

        // Perform simple translation along the direction from the center of muscle patch to bone
        translate_V(musVNew, t);
        int added = musVNew.rows();
        for (int i = 0; i < added; i++) {
            bb(cum_bnd + i) = Vmidx[i];
        }

        for (int i = 0; i < added; i++) {
            Bonedest.row(cum_bnd + i) = musVNew.row(i);
        }

        for (int i = 0; i < added; i++) {
            musclef.row(cum_bnd + i) = Vmpatch1.row(i);
        }
        cum_bnd = cum_bnd + added;

        // Using move_patch
        /*Eigen::MatrixXd musVNew;
        map_move_patch(Vmpatch1, Fmpatch1, Vpatch1, Fpatch1, musVNew);

        Eigen::VectorXi bnd;
        igl::boundary_loop(Fmpatch1,bnd);

        Eigen::MatrixXd bnd_uv;
        igl::map_vertices_to_circle(Vmpatch1,bnd,bnd_uv);

        Eigen::MatrixXd V_uv;
        igl::harmonic(Vmpatch1,Fmpatch1,bnd,bnd_uv,1,V_uv);

        std::cout << V_uv << std::endl;

        for (int i = 0; i < bnd.rows(); i++) {
            bb(cum_bnd + i) = Vmidx[bnd(i)];
        }

        for (int i = 0; i < bnd.rows(); i++) {
            Bonedest.row(cum_bnd + i) = musVNew.row(bnd(i));
        }

        for (int i = 0; i < bnd.rows(); i++) {
            musclef.row(cum_bnd + i) = Vmpatch1.row(bnd(i));
        }
        cum_bnd = cum_bnd + bnd.rows();*/
    }
    igl::ARAPData data;
    Eigen::VectorXi bbs = bb.segment(0, cum_bnd);
    igl::arap_precomputation(Vm, Fm, 3, bbs, data);
    //int ninterp = 10;
    int ninterp = 10;
    for (int ni = 1; ni <= ninterp; ni++) {
        double t = double (ni) / double (ninterp);
        //std::cout<< (1-t)*Fcc + t*Bcc <<std::endl;
        Eigen::MatrixXd currdest = (1-t)*musclef.block(0, 0, cum_bnd, 3) + t*Bonedest.block(0, 0, cum_bnd, 3);
        igl::arap_solve(currdest, data, Vm);
    }
//    double area1 = 0; double num1 = 0;
//    Eigen::MatrixXi Fpatch1(1, 3);
//    int Vidx[3*selected_faces1.size()];
//    std::set<int> s1;
//    int numdistinct = 0;
//    Eigen::RowVector3d center1 = Eigen::RowVector3d::Zero();
//
//    for (auto it = selected_faces1.begin(); it != selected_faces1.end(); ++it) {
//        // currently only one face is selected in face 1
//        ++it;
//        if (it == selected_faces1.end()) {
//            break;
//        }
//        --it;
//        // delete the above portion when selected_faces1 and selected_faces are correctly implemented
//
//        num1 += 1;
//        int curr = std::distance(selected_faces1.begin(), it);
//        Eigen::RowVectorXi triangle = F.row(*it);
//        Eigen::RowVector3d triangleCenter = (V.row(triangle(0)) + V.row(triangle(1)) + V.row(triangle(2))) / 3.0;
//        center1 += triangleCenter;
//
//        double a = (V.row(triangle(0)) - V.row(triangle(1))).norm();
//        double b = (V.row(triangle(1)) - V.row(triangle(2))).norm();
//        double c = (V.row(triangle(2)) - V.row(triangle(0))).norm();
//        double s = (a+b+c) / 2.;
//        area1 += pow(s*(s-a)*(s-b)*(s-c), 0.5);
//
//        // Construct smaller mesh for the patch on bone
//        for (int j = 0; j < 3; j++) {
//            int fj = triangle(j);
//            auto result = s1.insert(fj);
//            if ( result.second ) {
//                Vidx[numdistinct] = fj;
//                Fpatch1(curr, j) = numdistinct;
//                numdistinct++;
//            } else {
//                for (int k = 0; k < numdistinct; k++) {
//                    if (Vidx[k] == fj) {
//                        Fpatch1(curr, j) = k;
//                    }
//                }
//            }
//        }
//
//    }
//    center1 = 1./num1 * center1;
//    Eigen::MatrixXd Vpatch1(numdistinct, 3);
//    for (int i = 0; i < numdistinct; i++) {
//        Vpatch1.row(i) = V.row(Vidx[i]);
//    }
//
//    std::cout << "num1 " << num1 <<std::endl;
//    std::cout << "center1 " << center1 <<std::endl;
//
//
//
//    // Determine the total number of faces on the muscle to be deformed towards bone
//    int sourcenum1 = ceil(area1 / avg_area);
////    sourcenum1 = 100;
//    std::cout << "sourcenum1 " << sourcenum1 << std::endl;
//    Eigen::MatrixXd lst_dist(Fm.rows(), 1);
//    triangle_hunt_lst(center1, Vm, Fm, lst_dist);
//    Eigen::MatrixXd sY, sIX;
//    igl::sort(lst_dist, 1, true, sY, sIX);
//
//    Eigen::MatrixXd closest1 = sIX.block(0, 0, sourcenum1, 1);
//    Eigen::MatrixXi Fmpatch1(closest1.rows(), 3);
//    int Vmidx[3*closest1.rows()];
//    numdistinct = 0;
//
//    // Get the total number of distinct vertices
//    std::set<int> sm1;
//    for (int i = 0; i < closest1.rows(); i++) {
//        for (int j = 0; j < 3; j++) {
//            int fj = Fm(closest1(i), j);
//            auto result = sm1.insert(fj);
//            if ( result.second ) {
//                Vmidx[numdistinct] = fj;
//                Fmpatch1(i, j) = numdistinct;
//                numdistinct++;
//            } else {
//                for (int k = 0; k < numdistinct; k++) {
//                    if (Vmidx[k] == fj) {
//                        Fmpatch1(i, j) = k;
//                    }
//                }
//            }
//        }
//    }
//    std::cout << "loop " << std::endl;
//    Eigen::MatrixXd Vmpatch1(numdistinct, 3);
//    for (int i = 0; i < numdistinct; i++) {
//        Vmpatch1.row(i) = Vm.row(Vmidx[i]);
//    }
//    std::cout << "loop2 " << std::endl;
//
////    for (int i = 0; i < closest1.rows(); i++) {
////        Fpatch1.row(i) = Fm.row(closest1(i));
////    }
//
//    Eigen::Vector3d t = Eigen::Vector3d::Ones() * 1;
////    translate_V(Vpatch1, t);
//    std::cout << "loop 3" << std::endl;
//
//    Eigen::MatrixXd musVNew;
//    map_move_patch(Vmpatch1, Fmpatch1, Vpatch1, Fpatch1, musVNew);
//
//    std::cout << "loop 4" << std::endl;
//    Eigen::VectorXi bnd;
//    igl::boundary_loop(Fmpatch1,bnd);
//    std::cout << "loop 5" << std::endl;
//
////    VV.push_back(Vmpatch1);
////    FF.push_back(Fmpatch1);
////    VV.push_back(Vpatch1);
////    FF.push_back(Fpatch1);
//
//    Eigen::MatrixXd bnd_uv;
//    igl::map_vertices_to_circle(Vmpatch1,bnd,bnd_uv);
//
//    Eigen::MatrixXd V_uv;
//    igl::harmonic(Vmpatch1,Fmpatch1,bnd,bnd_uv,1,V_uv);
//
//    std::cout << V_uv << std::endl;
//
//    Eigen::VectorXi bb = Eigen::VectorXi(bnd.rows());
//    for (int i = 0; i < bnd.rows(); i++) {
//        bb(i) = Vmidx[bnd(i)];
//    }
//    Eigen::MatrixXd Bonedest(bnd.rows(), 3);
//    for (int i = 0; i < bnd.rows(); i++) {
//        Bonedest.row(i) = musVNew.row(bnd(i));
//    }
//    Eigen::MatrixXd musclef(bnd.rows(), 3);
//    for (int i = 0; i < bnd.rows(); i++) {
//        musclef.row(i) = Vmpatch1.row(bnd(i));
//    }
//    std::cout << "end " << std::endl;


//    igl::ARAPData data;
//    igl::arap_precomputation(Vm, Fm, 3, bb, data);
//    //int ninterp = 10;
//    int ninterp = 10;
//    for (int ni = 1; ni <= ninterp; ni++) {
//        double t = double (ni) / double (ninterp);
//        //std::cout<< (1-t)*Fcc + t*Bcc <<std::endl;
//        Eigen::MatrixXd currdest = (1-t)*musclef + t*Bonedest;
//        igl::arap_solve(currdest, data, Vm);
//    }
//
//
//    VV.push_back(Vm);
//    FF.push_back(Fm);
}
void construct_mesh(int * Vidx, const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    Eigen::MatrixXi & curr_patch,
                    Eigen::MatrixXd & Vpatch,
                    Eigen::MatrixXi & Fpatch) {
    Fpatch.resize(curr_patch.rows(), 3);
    std::set<int> s1;
    int numdistinct = 0;
    // Get the total number of distinct vertices
    for (int i = 0; i < curr_patch.rows(); i++) {
        Eigen::RowVector3i ci = curr_patch.row(i);
        for (int j = 0; j < 3; j++) {
            int fj = ci(j);
            auto result = s1.insert(fj);
            if ( result.second ) {
                Vidx[numdistinct] = fj;
                Fpatch(i, j) = numdistinct;
                numdistinct++;
            } else {
                for (int k = 0; k < numdistinct; k++) {
                    if (Vidx[k] == fj) {
                        Fpatch(i, j) = k;
                    }
                }
            }
        }
    }
    std::cout << "Fpatch finished" <<std::endl;
    Vpatch.resize(numdistinct, 3);
    for (int i = 0; i < numdistinct; i++) {
        std::cout << Vidx[i] << std::endl;
        Vpatch.row(i) = V.row(Vidx[i]);
    }
    std::cout << "Vpatch finished" <<std::endl;
}

void get_closest_patch(const Eigen::MatrixXd & V,
                        const Eigen::MatrixXi & F,
                        int sourcenum1,
                        Eigen::RowVector3d center1,
                        Eigen::MatrixXi & closest1) {
    Eigen::MatrixXd lst_dist(F.rows(), 1);
    triangle_hunt_lst(center1, V, F, lst_dist);
    Eigen::MatrixXd sY;
    Eigen::MatrixXi sIX;
    igl::sort(lst_dist, 1, true, sY, sIX);

    closest1 = sIX.block(0, 0, sourcenum1, 1);
}

void interpolate_ellipse(
        int npoints,
        struct elipse_param ep,
        struct elipse_param epm,
        Eigen::MatrixXd & curve,
        Eigen::MatrixXd & normal,
        Eigen::MatrixXd & long_dir,
        Eigen::MatrixXd & short_dir,
        Eigen::VectorXd & long_axis,
        Eigen::VectorXd & short_axis) {
    curve.resize(npoints, 3);
    normal.resize(npoints, 3);
    long_dir.resize(npoints, 3);
    short_dir.resize(npoints, 3);
    long_axis.resize(npoints);
    short_axis.resize(npoints);

    // Try to fix a little bit
    Eigen::Vector3d dir = epm.center - ep.center;
    epm.center += 0.1 * dir;
    ep.center += - 0.05 * dir;

    for (int i = 0; i < npoints; i++) {
        double t = (double) i/ (double) npoints;
        curve.row(i) = t * epm.center + (1-t) * ep.center;
        normal.row(i) = t * epm.normal.transpose() + (1-t) * ep.normal.transpose();
        long_dir.row(i) = t * epm.long_dir.transpose() + (1-t) * ep.long_dir.transpose();
        short_dir.row(i) = t * epm.short_dir.transpose() + (1-t) * ep.short_dir.transpose();
        long_axis(i) = (pow(t, 1) * epm.long_axis + (1-pow(t, 1)) * ep.long_axis) * (1.2*pow(t-0.5, 2) + 0.7 );
        short_axis(i) = (pow(t, 1) * epm.short_axis + (1-pow(t, 1)) * ep.short_axis ) * (1.2*pow(t-0.5, 2) + 0.7 ) ;
    }
}

void generate_tendon(const Eigen::MatrixXd & points,
                     const int num_points,
                     const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     const std::vector<Eigen::MatrixXi> & selected_faces,
                     Eigen::MatrixXd & Vm,
                     Eigen::MatrixXi & Fm,
                     Eigen::MatrixXd & Vt,
                     Eigen::MatrixXi & Ft) {
    Eigen::MatrixXd p  = points.block(points.rows() - num_points, 0, num_points, 3);
    p.colwise().reverse();
    // Your code here to populate V and F
    Eigen::MatrixXd Bc, Nc, pV, pN, Vm_before_smooth;
//    bezier(p1, p2, p3, 50, Bc, Nc);
//    bezier(p, num_points-1, 50, Bc, Nc);
    CatmullRomChain(p, 50, Bc, Nc);
//            bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, Bc, Nc);
    volume_along_curve(Bc, Nc, pV, pN);
    Eigen::MatrixXd All(pV.rows(), 6);
    All << pV, pN;
    poisson_surface_reconstruction(pV, pN, Vm, Fm);
    // Smooth before doing any of the face-matching operations
//    igl::per_vertex_attribute_smoothing(Vm_before_smooth, Fm, Vm);

    // Calculate average surface area of muscle mesh
    Eigen::MatrixXd dblA;
    igl::doublearea(Vm, Fm, dblA);
    double muscle_totalarea = 0;
    for (int i = 0; i < dblA.rows(); i++) {
        muscle_totalarea += dblA(i);
    }
    double avg_area = muscle_totalarea / dblA.rows();


    // Calculate area and center of attachment site 1 on bone
    // Also get the total number of distinct vertices in order
    // to construct a small mesh for the patch on bone.
    for (auto it = selected_faces.begin(); it != selected_faces.end(); ++it) {
        Eigen::MatrixXi curr_patch = *it;
        double area1 = 0; double num1 = 0;
        Eigen::MatrixXi Fpatch1;
        Eigen::MatrixXd Vpatch1;
        int Vidx[3*curr_patch.rows()];
        construct_mesh(Vidx, V, F, curr_patch, Vpatch1, Fpatch1);

        std::cout<<"construct mesh " <<std::endl;
        Eigen::RowVector3d center1 = Eigen::RowVector3d::Zero();
        int curr = std::distance(selected_faces.begin(), it);
//        if (curr == 1) {
//            break;
//        }
        for (int i = 0; i < curr_patch.rows(); i++) {
            num1 += 1;
            Eigen::RowVectorXi triangle = curr_patch.row(i);
            Eigen::RowVector3d triangleCenter = (V.row(triangle(0)) + V.row(triangle(1)) + V.row(triangle(2))) / 3.0;
            center1 += triangleCenter;

            double a = (V.row(triangle(0)) - V.row(triangle(1))).norm();
            double b = (V.row(triangle(1)) - V.row(triangle(2))).norm();
            double c = (V.row(triangle(2)) - V.row(triangle(0))).norm();
            double s = (a+b+c) / 2.;
            area1 += pow(s*(s-a)*(s-b)*(s-c), 0.5);
        }
        center1 = 1./num1 * center1;

        // Determine the total number of faces on the muscle to be connected and
        // generate a muscle patch
        int sourcenum1 = ceil(area1 / avg_area);
        std::cout << sourcenum1 <<std::endl;
        Eigen::MatrixXi closest1;
        get_closest_patch(Vm, Fm, sourcenum1, center1, closest1);
        std::cout<<"get closest patch m " <<std::endl;
        Eigen::MatrixXi Fmpatch1;
        Eigen::MatrixXd Vmpatch1;
        int Vmidx[Fm.rows()*3];
        std::cout<<closest1<<std::endl;
        Eigen::MatrixXi closest_patch(closest1.rows(), 3);
        // igl slice
        for (int i = 0; i < closest1.rows(); i++) {
            closest_patch.row(i) = Fm.row(closest1(i));
        }
        construct_mesh(Vmidx, Vm, Fm, closest_patch, Vmpatch1, Fmpatch1);

        std::cout<<"construct mesh  " <<std::endl;

        struct elipse_param ep = PCA_param(Vpatch1, Fpatch1);
        struct elipse_param epm = PCA_param(Vmpatch1, Fmpatch1);

        // Linearly interpolate the parameters of ellipse for poisson reconstruction
        Eigen::MatrixXd tendonV, tendonN;
        Eigen::MatrixXd curve, normal, long_dir, short_dir;
        Eigen::VectorXd long_axis, short_axis;

        // Fill in the interpolated values and generate point cloud
        interpolate_ellipse(20, ep, epm, curve, normal, long_dir, short_dir, long_axis, short_axis);
        ellipse_along_curve(curve, normal, long_axis, long_dir, short_axis, short_dir, tendonV, tendonN);

        Eigen::MatrixXd Vt_before_smooth;
        poisson_surface_reconstruction(tendonV, tendonN, Vt_before_smooth, Ft);
        igl::per_vertex_attribute_smoothing(Vt_before_smooth, Ft, Vt);

        Eigen::Vector3d t = Eigen::Vector3d::Ones() * 1;
        translate_V(Vmpatch1, t);

        std::cout << "elipse ep print" << std::endl;
        std::cout << ep.long_axis << std::endl;
        std::cout << ep.long_dir << std::endl;
        std::cout << ep.short_axis << std::endl;
        std::cout << ep.short_dir << std::endl;
        std::cout << ep.center << std::endl;
        std::cout << ep.normal << std::endl;

        std::cout << "elipse epm print" << std::endl;
        std::cout << epm.long_axis << std::endl;
        std::cout << epm.long_dir << std::endl;
        std::cout << epm.short_axis << std::endl;
        std::cout << epm.short_dir << std::endl;
        std::cout << epm.center << std::endl;
        std::cout << epm.normal << std::endl;
        Eigen::MatrixXd musVNew;


//        map_move_patch(Vmpatch1, Fmpatch1, Vpatch1, Fpatch1, musVNew);
//
//        Eigen::VectorXi bnd;
//        igl::boundary_loop(Fmpatch1,bnd);
//
//        Eigen::MatrixXd bnd_uv;
//        igl::map_vertices_to_circle(Vmpatch1,bnd,bnd_uv);
//
//        Eigen::MatrixXd V_uv;
//        igl::harmonic(Vmpatch1,Fmpatch1,bnd,bnd_uv,1,V_uv);
//
//        std::cout << V_uv << std::endl;
//
//        Eigen::VectorXi bb = Eigen::VectorXi(bnd.rows());
//        for (int i = 0; i < bnd.rows(); i++) {
//            bb(i) = Vmidx[bnd(i)];
//        }
//        Eigen::MatrixXd Bonedest(bnd.rows(), 3);
//        for (int i = 0; i < bnd.rows(); i++) {
//            Bonedest.row(i) = musVNew.row(bnd(i));
//        }
//        Eigen::MatrixXd musclef(bnd.rows(), 3);
//        for (int i = 0; i < bnd.rows(); i++) {
//            musclef.row(i) = Vmpatch1.row(bnd(i));
//        }
//        std::cout << "end " << std::endl;

    }
}
//    igl::ARAPData data;
//    igl::arap_precomputation(Vm, Fm, 3, bb, data);
//    //int ninterp = 10;
//    int ninterp = 10;
//    for (int ni = 1; ni <= ninterp; ni++) {
//        double t = double (ni) / double (ninterp);
//        //std::cout<< (1-t)*Fcc + t*Bcc <<std::endl;
//        Eigen::MatrixXd currdest = (1-t)*musclef + t*Bonedest;
//        igl::arap_solve(currdest, data, Vm);
//    }
