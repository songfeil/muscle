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

int main(int argc, char *argv[])
{
  std::vector<Eigen::MatrixXd> VV; // vector of vertex matrices for all meshes
  std::vector<Eigen::MatrixXi> FF; // vector of face matrices for all meshes
  // Load in a mesh
//  igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/knight.off", V, F);

//
//  Eigen::MatrixXi E = edges(F);
//  int Chi = euler_characteristic(F);
//  std::cout<<"Edge list E is "<<E.rows()<<"x"<<E.cols()<<std::endl;
//  std::cout<<"Euler Characteristic: "<<Chi<<std::endl;

  // OUR CODE HERE!
  //generate_bone(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 3, 5), V, F);

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

  // Mode consts
  enum Mode
  {
    NONE = 0,
    BONE = 1,
    MUSCLE = 2,
  } mode = NONE;

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
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        igl::combine(VV,FF,V,F);
        viewer.data().set_mesh(V,F);
      }
    }
    else if (mode == MUSCLE) {
      viewer.data().clear();
      viewer.data().set_points(muscle_points, blue);
      if (VV.size() > 0) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        igl::combine(VV,FF,V,F);
        viewer.data().set_mesh(V,F);
      }
    }
  };

  viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    if (mode > 0) {
      // Ray cast to x-y plane
      last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);
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
      if (mode == BONE) {
        // Add intersection to list of bone points
        bone_points.conservativeResize(bone_points.rows() + 1, 3);
        bone_points.row(bone_points.rows() - 1) = intersection;
        n_bp++;
        // For every 2 points, generate a bone
        if (n_bp == 2){
          Eigen::MatrixXd V;
          Eigen::MatrixXi F;
          Eigen::Vector3d p1 = bone_points.row(bone_points.rows() - 1);
          Eigen::Vector3d p2 = bone_points.row(bone_points.rows() - 2);
          generate_bone(p1, p2, V, F);
          // Add to list of meshes
          VV.push_back(V);
          FF.push_back(F);
          n_bp = 0;
        }
      }
        else if (mode == MUSCLE) {
          // Add intersection to list of muscle points
          muscle_points.conservativeResize(muscle_points.rows() + 1, 3);
          muscle_points.row(muscle_points.rows() - 1) = intersection;
          n_cp++;
          // For every 3 points, generate a muscle!
          if (n_cp == 3) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Eigen::Vector3d p1 = muscle_points.row(muscle_points.rows() - 1);
            Eigen::Vector3d p2 = muscle_points.row(muscle_points.rows() - 2);
            Eigen::Vector3d p3 = muscle_points.row(muscle_points.rows() - 3);
            // Your code here to populate V and F
            Eigen::MatrixXd Bc, Nc, pV, pN;
            bezier(p1, p2, p3, 4, Bc, Nc);
//            bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, Bc, Nc);
            volume_along_curve(Bc, Nc, pV, pN);
            std::cout << pV << std::endl;
            poisson_surface_reconstruction(pV, pN, V, F);
//            deform(p1, p2, p3, V, F);
            VV.push_back(V);
            FF.push_back(F);
            std::cout << "generate muscle" << std::endl;
            n_cp = 0;
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
    }
    update();
    return true;
  };
    

  // Launch a viewer instance
  viewer.launch();
  return 0;
}

