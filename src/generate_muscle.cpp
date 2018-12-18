//
// Created by 宋飞龙 on 2018-11-11.
//

#include <iostream>
#include <igl/arap.h>
#include "generate_muscle.h"
#include "control_curve.h"
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
#include <PCA_ellipse.h>
#include <igl/per_vertex_attribute_smoothing.h>
#include <igl/rotation_matrix_from_directions.h>

// Translate vertices V by vector t
void translate_V(Eigen::MatrixXd & V, const Eigen::Vector3d & t) {
    for (int i = 0; i < V.rows(); i++) {
        V.row(i) += t.transpose();
    }
}

// Construct patch mesh originally coming from a larger mesh,
// with all indices into V, F converted to local indexing
// V, F: larger mesh where the patch mesh comes from
// curr_patch: indices into rows of F that indicate faces included
//             in the small patch mesh
// Vpatch, Fpatch: output patch mesh with local indexing
// Vidx: list that maps the local index of vertices to the original global index
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
    Vpatch.resize(numdistinct, 3);
    for (int i = 0; i < numdistinct; i++) {
        Vpatch.row(i) = V.row(Vidx[i]);
    }
}

// Find top k closest patches on the mesh (V,F) from center
void get_closest_patch(const Eigen::MatrixXd & V,
                        const Eigen::MatrixXi & F,
                        int k,
                        Eigen::RowVector3d center,
                        Eigen::MatrixXi & closest1) {
    Eigen::MatrixXd lst_dist(F.rows(), 1);
    triangle_hunt_lst(center, V, F, lst_dist);
    Eigen::MatrixXd sY;
    Eigen::MatrixXi sIX;
    igl::sort(lst_dist, 1, true, sY, sIX);

    closest1 = sIX.block(0, 0, k, 1);
}

// Interpolate the ellipse parameters
void interpolate_ellipse(
        int npoints,
        struct ellipse_param ep,
        struct ellipse_param epm,
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

void generate_muscle(const Eigen::MatrixXd & points,
                     const int num_points,
                     Eigen::MatrixXd & Vm,
                     Eigen::MatrixXi & Fm) {

    Eigen::MatrixXd Vm_before_smooth;
    //Organize our input point data
    Eigen::MatrixXd p  = points.block(points.rows() - num_points, 0, num_points, 3);
    p.colwise().reverse();
    // Fit spline to input points
    Eigen::MatrixXd Bc, Nc, pV, pN;
    CatmullRomChain(p, 50, Bc, Nc);
    // Use spline to generate point cloud
    volume_along_curve(Bc, Nc, pV, pN);
    // Reconstruct surface mesh from point cloud, marching cubes
    poisson_surface_reconstruction(pV, pN, Vm_before_smooth, Fm);
    // Smooth (marching cubes is ugly)
    igl::per_vertex_attribute_smoothing(Vm_before_smooth, Fm, Vm);
}

void attach_muscle( const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    Eigen::MatrixXd & Vm,
                    Eigen::MatrixXi & Fm,
                    std::set<int> & attached_vids) {

    // Matrices containing info on the vertices we are moving
    Eigen::VectorXi bb = Eigen::VectorXi(selected_faces.size() * 3);
    Eigen::MatrixXd Bonedest(selected_faces.size()*3, 3);
    Eigen::MatrixXd musclef(selected_faces.size()*3, 3);

    // Matrices containing info on already fixed vertices
    Eigen::VectorXi bb_fixed = Eigen::VectorXi(attached_vids.size());
    Eigen::MatrixXd fixed_verts(attached_vids.size(), 3);

    // Ensure to fix the already attached vertices before adding new attached vertices to the set
    for (auto it = attached_vids.begin(); it != attached_vids.end(); ++it){
        int vid = *it;
        int curr = std::distance(attached_vids.begin(), it);
        bb_fixed(curr) = vid;
        fixed_verts.row(curr) = Vm.row(vid);
    }

    for (auto it = selected_faces.begin(); it != selected_faces.end(); ++it) {
        int curr = std::distance(selected_faces.begin(), it);
        // Our current selected triangle on the bone
        Eigen::RowVectorXi triangle = (*it).row(0);
        Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
        P.row(0) = V.row(triangle(0));
        P.row(1) = V.row(triangle(1));
        P.row(2) = V.row(triangle(2));
        // Find the closest triangle on the muscle
        int fi = triangle_hunt(P, Vm, Fm);
        Eigen::RowVector3i f = Fm.row(fi);
        attached_vids.insert(f(0));
        attached_vids.insert(f(1));
        attached_vids.insert(f(2));
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
    }
    
    Eigen::VectorXi bb_all = Eigen::VectorXi(bb.rows() + bb_fixed.rows());
    if (bb_fixed.rows() > 0){
        bb_all << bb, bb_fixed;
    }
    else {
        bb_all = bb;
    }
    igl::ARAPData data;
    igl::arap_precomputation(Vm, Fm, 3, bb_all, data);
    int ninterp = 10;
    for (int ni = 1; ni <= ninterp; ni++) {
        double t = double (ni) / double (ninterp);
        Eigen::MatrixXd currdest = (1-t)*musclef + t*Bonedest;
        Eigen::MatrixXd all_verts(currdest.rows() + fixed_verts.rows(), 3);
        if (fixed_verts.rows() > 0) {
            all_verts << currdest, fixed_verts;
        }
        else {
            all_verts = currdest;
        }
        igl::arap_solve(all_verts, data, Vm);
    }
}

// Convert cartesian coordinate to polar coordinate
void card2pol(Eigen::MatrixXd & PCAproj, Eigen::MatrixXd & PCApol) {
    for (int i = 0; i < PCApol.rows(); i++) {
        double r = PCAproj.row(i).norm();
        PCApol(i, 1) = r;
        if (r > 0) {
            if (PCAproj(i, 1) > 0) {
                PCApol(i, 0) = acos(PCAproj(i, 0) / r);
            } else {
                PCApol(i, 0) = 2*M_PI - acos(PCAproj(i, 0) / r);
            }
        } else {
            PCApol(i, 0) = 0;
        }
    }
}

void attach_muscle_multiface( const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    Eigen::MatrixXd & Vm,
                    Eigen::MatrixXi & Fm,
                    std::set<int> & attached_vids) {

    int largenum = 10000;
    Eigen::VectorXi bb = Eigen::VectorXi(largenum);
    Eigen::MatrixXd Bonedest(largenum, 3);
    Eigen::MatrixXd musclef(largenum, 3);

    // Matrices containing info on already fixed vertices
    Eigen::VectorXi bb_fixed = Eigen::VectorXi(attached_vids.size());
    Eigen::MatrixXd fixed_verts(attached_vids.size(), 3);

    // Ensure to fix the already attached vertices before adding new attached vertices to the set
    for (auto it = attached_vids.begin(); it != attached_vids.end(); ++it){
        int vid = *it;
        int curr = std::distance(attached_vids.begin(), it);
        bb_fixed(curr) = vid;
        fixed_verts.row(curr) = Vm.row(vid);
    }

    // Calculate average surface area of muscle mesh
    Eigen::MatrixXd dblA;
    igl::doublearea(Vm, Fm, dblA);
    double muscle_totalarea = 0;
    for (int i = 0; i < dblA.rows(); i++) {
        muscle_totalarea += dblA(i);
    }
    double avg_area = muscle_totalarea / dblA.rows();

    int curr_num = 0;

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

        Eigen::RowVector3d center1 = Eigen::RowVector3d::Zero();
        int curr = std::distance(selected_faces.begin(), it);
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
        Eigen::MatrixXi closest1;
        get_closest_patch(Vm, Fm, sourcenum1, center1, closest1);
        Eigen::MatrixXi Fmpatch1;
        Eigen::MatrixXd Vmpatch1;
        int Vmidx[Fm.rows()*3];
        Eigen::MatrixXi closest_patch(closest1.rows(), 3);
        // igl slice
        for (int i = 0; i < closest1.rows(); i++) {
            closest_patch.row(i) = Fm.row(closest1(i));
        }
        construct_mesh(Vmidx, Vm, Fm, closest_patch, Vmpatch1, Fmpatch1);

        // Adding the vertices of the new muscle attachment points to our list
        // So that future deformations ensure these vertices remain fixed
        for (int i = 0; i < closest_patch.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                attached_vids.insert(closest_patch(i, j));
            }
        }

        struct ellipse_param ep = PCA_param(Vpatch1, Fpatch1);
        struct ellipse_param epm = PCA_param(Vmpatch1, Fmpatch1);

        // Project to PCA basis
        ep.long_dir.normalize();
        ep.short_dir.normalize();
        ep.normal.normalize();
        epm.long_dir.normalize();
        epm.short_dir.normalize();
        epm.normal.normalize();

        // Calculate average normal of the patches and align the last principal
        // direction to the average normal direction
        Eigen::MatrixXd bonN, musN;
        igl::per_face_normals(Vpatch1, Fpatch1, bonN);
        igl::per_face_normals(Vmpatch1, Fmpatch1, musN);

        // Average per face normal
        Eigen::VectorXd avgbonN, avgmusN;
        avgbonN = bonN.colwise().mean();
        avgmusN = musN.colwise().mean();

        Eigen::Vector3d bon_vecnormal = ep.normal.col(0);
        Eigen::Vector3d mus_vecnormal = epm.normal.col(0);

        // Flip the normal direction if needed
        if (bon_vecnormal.dot(avgbonN) < 0) {
            bon_vecnormal = - bon_vecnormal;
        }
        if (mus_vecnormal.dot(avgmusN) > 0) {
            mus_vecnormal = - mus_vecnormal;
        }

        // Rotate the bone patch by aligning the patch normal to that of muscle patch
        Eigen::Matrix3d rotation_mat = igl::rotation_matrix_from_directions(bon_vecnormal, mus_vecnormal);

        // Translate the bone patch center to the muscle patch center
        Eigen::MatrixXd Vpatch1_centered(Vpatch1.rows(), 3);
        for (int i = 0 ; i < Vpatch1.rows(); i++) {
            Vpatch1_centered.row(i) = Vpatch1.row(i);
        }
        translate_V(Vpatch1_centered, -ep.center);
        Eigen::MatrixXd Vmpatch1_centered(Vmpatch1.rows(), 3);
        for (int i = 0 ; i < Vmpatch1.rows(); i++) {
            Vmpatch1_centered.row(i) = Vmpatch1.row(i);
        }
        translate_V(Vmpatch1_centered, -epm.center);
        Eigen::MatrixXd Vpatch1_trans = Vpatch1_centered * rotation_mat.transpose();

        // Get right-hand sided standard coordinate system by fixing the directions of the PCA basis
        if (mus_vecnormal.dot(epm.long_dir.cross(epm.short_dir)) < 0) {
            epm.short_dir = -epm.short_dir;
        }

        // Project the patch points onto the plane spanned by the first two principal directions
        Eigen::MatrixXd bonPCAproj(Vpatch1.rows(), 2);
        Eigen::MatrixXd musPCAproj(Vmpatch1.rows(), 2);
        bonPCAproj.col(0) = Vpatch1_trans * epm.long_dir;
        bonPCAproj.col(1) = Vpatch1_trans * epm.short_dir;
        musPCAproj.col(0) = Vmpatch1_centered * epm.long_dir;
        musPCAproj.col(1) = Vmpatch1_centered * epm.short_dir;

        // Transform to polar coordinates
        Eigen::MatrixXd bonPCApol(Vpatch1.rows(), 2);
        Eigen::MatrixXd musPCApol(Vmpatch1.rows(), 2);
        card2pol(bonPCAproj, bonPCApol);
        card2pol(musPCAproj, musPCApol);

        // Get the parametrization of the patch points
        Eigen::VectorXi bonBND;
        igl::boundary_loop(Fpatch1,bonBND);
        Eigen::MatrixXd bonbnd_uv;
        igl::map_vertices_to_circle(Vmpatch1,bonBND,bonbnd_uv);
        Eigen::MatrixXd bonV_uv;
        igl::harmonic(Vpatch1,Fpatch1,bonBND,bonbnd_uv,1,bonV_uv);

        // Flipping the connection matrix so that all normals are flipped
        Eigen::MatrixXi Fmpatch1_flip(Fmpatch1.rows(), 3);
        for (int i = 0; i < Fmpatch1.rows(); i++) {
            Fmpatch1_flip(i, 2) = Fmpatch1(i, 0);
            Fmpatch1_flip(i, 1) = Fmpatch1(i, 1);
            Fmpatch1_flip(i, 0) = Fmpatch1(i, 2);
        }

        Eigen::VectorXi musBND;
        igl::boundary_loop(Fmpatch1_flip,musBND);
        Eigen::MatrixXd musbnd_uv;
        igl::map_vertices_to_circle(Vmpatch1,musBND,musbnd_uv);
        Eigen::MatrixXd musV_uv;
        igl::harmonic(Vmpatch1,Fmpatch1_flip,musBND,musbnd_uv,1,musV_uv);

        // Transform to polar coordinate system
        Eigen::MatrixXd bonV_uvpol(Vpatch1.rows(), 2);
        Eigen::MatrixXd musV_uvpol(Vmpatch1.rows(), 2);
        card2pol(bonV_uv, bonV_uvpol);
        card2pol(musV_uv, musV_uvpol);

        // Find the optimal rotation angle to align project on PCA and parametrization
        // by minimizing (thetaproj - (thetaparam + t)) ^ 2. Note thetaproj - thetaparam
        // should be the angle between the projection point and parametrization point
        Eigen::MatrixXd bonmeanuvpol = bonV_uvpol.colwise().mean();
        Eigen::MatrixXd musmeanuvpol = musV_uvpol.colwise().mean();
        double t_bon_uv = bonmeanuvpol(0);
        double t_mus_uv = musmeanuvpol(0);

        Eigen::MatrixXd bonmeanPCApol = bonPCApol.colwise().mean();
        Eigen::MatrixXd musmeanPCApol = musPCApol.colwise().mean();

        double t_bon_PCA = 0;
        double t_mus_PCA = 0;

        double t_bon = 0;
        double t_mus = 0;
        for (int i = 0; i < bonV_uvpol.rows(); i++) {
            double thetaproj = bonPCApol(i, 0);
            double thetaparam = bonV_uvpol(i, 0);
            double diff = thetaproj - thetaparam;
            if (diff < 0) {
                // Thetaproj and thetaparam crosses angle zero
                diff += 2*M_PI;
            }
            if (diff > M_PI) {
                diff = diff - 2*M_PI;
            }
            t_bon += diff;
        }
        // Minimizer for bone patch
        t_bon = t_bon / bonV_uvpol.rows();
        for (int i = 0; i < musV_uvpol.rows(); i++) {
            double thetaproj = musPCApol(i, 0);
            double thetaparam = musV_uvpol(i, 0);
            double diff = thetaproj - thetaparam;
            if (diff < 0) {
                // Thetaproj and thetaparam crosses angle zero
                diff += 2*M_PI;
            }
            if (diff > M_PI) {
                diff = diff - 2*M_PI;
            }
            t_mus += diff;
        }
        // Minimizer for bone patch
        t_mus = t_mus / musV_uvpol.rows();

        for (int i = 0; i < bonV_uvpol.rows(); i++) {
            bonV_uvpol(i, 0) += t_bon;
            if (bonV_uvpol(i, 0) < 0) {
                bonV_uvpol(i, 0) += 2*M_PI;
            }
            if (bonV_uvpol(i, 0) > 2*M_PI) {
                bonV_uvpol(i, 0) -= 2*M_PI;
            }
        }
        for (int i = 0; i < musV_uvpol.rows(); i++) {
            musV_uvpol(i, 0) += t_mus;
            if (musV_uvpol(i, 0) < 0) {
                musV_uvpol(i, 0) += 2*M_PI;
            }

            if (musV_uvpol(i, 0) > 2*M_PI) {
                musV_uvpol(i, 0) -= 2*M_PI;
            }
        }

        // Align boundary points on the muscle patch to the boundary points on the bone patch
        Eigen::MatrixXd musVNew(musBND.rows(), 3);
        for (int i = 0; i < musBND.rows(); i++) {
            bool found = false;
            for (int j = 0; j < bonBND.rows(); j++) {
                double bontheta0 = bonV_uvpol(bonBND(j), 0);
                double bontheta1 = bonV_uvpol(bonBND((j+1)%bonBND.rows()), 0);
                double mustheta  = musV_uvpol(musBND(i), 0);
                bool in = false;
                if (bontheta1 < bontheta0) {
                    bontheta1 += 2*M_PI;
                    if (mustheta >= bontheta0) {
                        in = true;
                    } else {
                        mustheta += 2*M_PI;
                        if (mustheta <= bontheta1) {
                            in = true;
                        }
                    }
                } else {
                    if ( bontheta0 <= mustheta && bontheta1 > mustheta) {
                        in = true;
                    }
                }

                if ( in ) {
                    found = true;
                    double t = (mustheta - bontheta0) / (bontheta1 - bontheta0);
                    musVNew.row(i) = t * Vpatch1.row(bonBND((j+1)%bonBND.rows())) + (1-t) * Vpatch1.row(bonBND(j));
                    break;
                }
            }
        }

        // Record the destination points of the muscle patch for deformation
        for (int i = 0; i < musBND.rows(); i++) {
            bb(curr_num +i) = Vmidx[musBND(i)];
        }
        for (int i = 0; i < musBND.rows(); i++) {
            Bonedest.row(curr_num +i) = musVNew.row(i);
        }
        for (int i = 0; i < musBND.rows(); i++) {
            musclef.row(curr_num +i) = Vmpatch1.row(musBND(i));
        }
        curr_num += musBND.rows();

    }

    Eigen::VectorXi bb_all = Eigen::VectorXi(bb.rows() + bb_fixed.rows());
    if (bb_fixed.rows() > 0){
        bb_all << bb.segment(0, curr_num), bb_fixed;
    }
    else {
        bb_all = bb.segment(0, curr_num);
    }
    igl::ARAPData data;
    igl::arap_precomputation(Vm, Fm, 3, bb_all, data);
    int ninterp = 10;
    for (int ni = 1; ni <= ninterp; ni++) {
        double t = double (ni) / double (ninterp);
        Eigen::MatrixXd currdest = (1-t)*musclef.block(0, 0, curr_num ,3) + t*Bonedest.block(0, 0, curr_num ,3);
        Eigen::MatrixXd all_verts(currdest.rows() + fixed_verts.rows(), 3);
        if (fixed_verts.rows() > 0) {
            all_verts << currdest, fixed_verts;
        }
        else {
            all_verts = currdest;
        }
        igl::arap_solve(all_verts, data, Vm);
    }
}

void attach_tendon(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<Eigen::MatrixXi> & selected_faces,
                    const Eigen::MatrixXd & Vm,
                    const Eigen::MatrixXi & Fm,
                    std::vector<Eigen::MatrixXd> & VV,
                    std::vector<Eigen::MatrixXi> & FF,
                    std::set<int> & attached_vids) {

    // Calculate average surface area of muscle mesh
    Eigen::MatrixXd Vt;
    Eigen::MatrixXi Ft;

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

        Eigen::RowVector3d center1 = Eigen::RowVector3d::Zero();
        int curr = std::distance(selected_faces.begin(), it);
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
        Eigen::MatrixXi closest1;
        get_closest_patch(Vm, Fm, sourcenum1, center1, closest1);
        Eigen::MatrixXi Fmpatch1; 
        Eigen::MatrixXd Vmpatch1;
        int Vmidx[Fm.rows()*3];
        Eigen::MatrixXi closest_patch(closest1.rows(), 3);
        // igl slice
        for (int i = 0; i < closest1.rows(); i++) {
            closest_patch.row(i) = Fm.row(closest1(i));
        }
        construct_mesh(Vmidx, Vm, Fm, closest_patch, Vmpatch1, Fmpatch1);

        // Adding the vertices of the new muscle attachment points to our list
        // So that future deformations ensure these vertices remain fixed
        for (int i = 0; i < closest_patch.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                attached_vids.insert(closest_patch(i, j));
            }
        }

        struct ellipse_param ep = PCA_param(Vpatch1, Fpatch1);
        struct ellipse_param epm = PCA_param(Vmpatch1, Fmpatch1);

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

        VV.push_back(Vt);
        FF.push_back(Ft);
    }
}
