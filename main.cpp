#include "edges.h"
#include "euler_characteristic.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include "generate_bone.h"
#include "bezier.h"
#include "ray_intersect_plane.h"
#include <igl/unproject.h>
#include <igl/unproject_ray.h>
#include <igl/combine.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/Geometry>
#include "cylinder.h"
#include "deform.h"
#include "pick_constrain_point.h"
#include <vector>
#include <iostream>
#include "generate_muscle.h"
#include "gaussian.h"
#include "volume_along_curve.h"
#include "poisson_surface_reconstruction.h"
#include <Eigen/Sparse>
#include <igl/unproject_onto_mesh.h>
#include "triangle_hunt.h"
#include <igl/arap.h>

// Mode consts
enum Mode
{
  NONE = 0,
  BONE = 1,
  MUSCLE = 2,
  FACE_SELECT = 3,
};

std::vector<Eigen::MatrixXd> VV;
std::vector<Eigen::MatrixXi> FF;
Eigen::MatrixXd V;
Eigen::MatrixXi F;
std::set<int> selected_faces;
Eigen::MatrixXd face_colors;
Mode mode;

int main(int argc, char *argv[])
{
   // vector of vertex matrices for all meshes
   // vector of face matrices for all meshes
  // Load in a mesh
//  igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/knight.off", V, F);

//
//  Eigen::MatrixXi E = edges(F);
//  int Chi = euler_characteristic(F);
//  std::cout<<"Edge list E is "<<E.rows()<<"x"<<E.cols()<<std::endl;
//  std::cout<<"Euler Characteristic: "<<Chi<<std::endl;

  // OUR CODE HERE!
  //generate_bone(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 3, 5), V, F);
  mode = NONE;

  Eigen::VectorXd G;
  gaussian(10, 1, 3, 1, G);
  gaussian(10, 1, 3, 2, G);

  Eigen::MatrixXd B, N;
  bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, B, N);


//  std::cout << B << std::endl;
//   Eigen::MatrixXd CV;
//   Eigen::MatrixXi CF;
// //  cylinder(12, 22, CV, CF);
//   deform(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(5, 5, 0), Eigen::Vector3d(0, 5, 0), CV, CF);
//   VV.push_back(CV);
//   FF.push_back(CF);

  // Create a libigl Viewer object 
  igl::opengl::glfw::Viewer viewer;
  std::cout << "========= MUSCLE/BONE GENERATION =========" << std::endl;
  std::cout << " B,b \t Bone input mode \n M,m \t Muscle input mode \n[space]  No mode (can drag around scene again)" << std::endl;

  // Color consts
  const Eigen::RowVector3d orange(1.0,0.7,0.2);
  const Eigen::RowVector3d yellow(1.0,0.9,0.2);
  const Eigen::RowVector3d blue(0.2,0.3,0.8);
  const Eigen::RowVector3d green(0.2,0.6,0.3);

  // User input data
  Eigen::Vector3d last_mouse;
  Eigen::MatrixXd bone_points;
  Eigen::MatrixXd muscle_points;
  int n_bp = 0;
  int n_cp = 0;

  // Define update functions for viewer
  const auto & update = [&]()
  {
    if (mode == BONE) {
      viewer.data().clear();
      viewer.data().set_points(bone_points, orange);
      // Setting multiple meshes
      // http://www.alecjacobson.com/weblog/?p=4679
      if (VV.size() > 0) {
        igl::combine(VV,FF,V,F);
        face_colors.resize(F.rows(), 3);
        for (int i = 0; i < F.rows(); i++) {
          face_colors.row(i) = green;
        }
        viewer.data().set_mesh(V,F);
        viewer.data().set_colors(face_colors);
      }
    }
    else if (mode == MUSCLE) {
      viewer.data().clear();
      viewer.data().set_points(muscle_points, blue);
      if (VV.size() > 0) {
        igl::combine(VV,FF,V,F);
        face_colors.resize(F.rows(), 3);
        for (int i = 0; i < F.rows(); i++) {
          face_colors.row(i) = green;
        }
        viewer.data().set_mesh(V,F);
        viewer.data().set_colors(face_colors);
      }
    }
    else if ( mode == FACE_SELECT) {
      viewer.data().set_colors(face_colors);
    }
  };

  viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    if ( mode > 0) {
      // Ray cast to x-y plane
      last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);
      if ( mode == FACE_SELECT) {
        // Find closest point on mesh to mouse position
        int fid;
        Eigen::Vector3f bary;
        Eigen::RowVector2f mouse;
        mouse << last_mouse.head(2)(0), last_mouse.head(2)(1);
        if(igl::unproject_onto_mesh(
          mouse,
          viewer.core.view,
          viewer.core.proj, 
          viewer.core.viewport, 
          V, F, 
          fid, bary))
        {
          // selected face! color it
          selected_faces.insert(fid);
          face_colors.row(fid) = orange;
        }
      }
      else {
        Eigen::Vector3d source, dir, intersection;
      Eigen::Vector2d p;
      p << last_mouse.head(2)(0), last_mouse.head(2)(1);
      igl::unproject_ray(
        p,
        viewer.core.view,
        viewer.core.proj, 
        viewer.core.viewport, 
        source, dir);
      ray_intersect_plane(Eigen::Vector3d(0, 0, 0),
                          Eigen::Vector3d(0, 0, 1),
                          source,
                          dir,
                          intersection);
      if ( mode == BONE) {
        // Add intersection to list of bone points
        bone_points.conservativeResize(bone_points.rows() + 1, 3);
        bone_points.row(bone_points.rows() - 1) = intersection;
        // n_bp++;
        // // For every 2 points, generate a bone
        // if (n_bp == 2){
        //   Eigen::MatrixXd V;
        //   Eigen::MatrixXi F;
        //   Eigen::Vector3d p1 = bone_points.row(bone_points.rows() - 1);
        //   Eigen::Vector3d p2 = bone_points.row(bone_points.rows() - 2);
        //   generate_bone(p1, p2, V, F);
        //   // Add to list of meshes
        //   VV.push_back(V);
        //   FF.push_back(F);
        //   n_bp = 0;
          // if (VV.size() == 2) {
          //   Eigen::MatrixXd Vnew;
          //   Eigen::MatrixXi Fnew;
          //   igl::copyleft::cgal::mesh_boolean(VV.at(0),FF.at(0),VV.at(1),FF.at(1),FB,MESH_BOOLEAN_TYPE_MINUS,Vnew,Fnew);
          //   VV.clear();
          //   FF.clear();
          //   VV.push_back(Vnew);
          //   FF.push_back(Fnew);
          ///}
        //}
      }
        else if ( mode == MUSCLE) {
          // Add intersection to list of muscle points
          muscle_points.conservativeResize(muscle_points.rows() + 1, 3);
          muscle_points.row(muscle_points.rows() - 1) = intersection;
          n_cp++;
          // For every 3 points, generate a muscle!
          if (n_cp == 3) {
            Eigen::MatrixXd Vm;
            Eigen::MatrixXi Fm;
            Eigen::Vector3d p1 = muscle_points.row(muscle_points.rows() - 1);
            Eigen::Vector3d p2 = muscle_points.row(muscle_points.rows() - 2);
            Eigen::Vector3d p3 = muscle_points.row(muscle_points.rows() - 3);
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
            n_cp = 0;
          }
        }
      }
        update();
        return true;
      }
    return false;
  };

  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch(key) {
      case 'B':
      case 'b':
      {
         mode = BONE;
        break;
      }
      case 'M':
      case 'm':
      {
         mode = MUSCLE;
        break;
      }
      case ' ':
      {
         mode = NONE;
        break;
      }
      case 'f':
      {
         mode = FACE_SELECT;
        break;
      }
      case 'G':
      case 'g':
      {
        if ( mode == BONE) {
          generate_bones(bone_points, VV, FF);
        }
        break;
      }
    }
    update();
    return true;
  };

  // Launch a viewer instance
  viewer.launch();
  return 0;
}

