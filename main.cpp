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
#include <iostream>

int main(int argc, char *argv[])
{
  std::vector<Eigen::MatrixXd> VV;
  std::vector<Eigen::MatrixXi> FF;
  std::vector<Eigen::MatrixXd> BB;

  // Create a libigl Viewer object 
  igl::opengl::glfw::Viewer viewer;

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
    CURVE = 2,
  } mode = CURVE;

  // User input data
  Eigen::Vector3d last_mouse;
  Eigen::MatrixXd bone_points;
  Eigen::MatrixXd curve_points;
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
    else if (mode == CURVE) {
      viewer.data().set_points(curve_points, blue);
    }
  };

  viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    if (mode > 0) {
      std::cout << "mode: " << mode << std::endl;
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
        else if (mode == CURVE) {
          // Add intersection to list of curve points
          std::cout << "curve?" << std::endl;
          curve_points.conservativeResize(curve_points.rows() + 1, 3);
          curve_points.row(curve_points.rows() - 1) = intersection;
          std::cout << curve_points << std::endl;
          n_cp++;
          // For every 3 points, generate a curve
          if (n_cp == 3) {
            Eigen::Vector3d p1 = curve_points.row(curve_points.rows() - 1);
            Eigen::Vector3d p2 = curve_points.row(curve_points.rows() - 2);
            Eigen::Vector3d p3 = curve_points.row(curve_points.rows() - 3);
            Eigen::MatrixXd B;
            bezier(p1, p2, p3, 10, B);
            BB.push_back(B);
            n_cp = 0;
            std::cout << "B \n" << B << std::endl;
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
      case 'b':
      {
        mode = BONE;
        break;
      }
      case 'c':
      {
        mode = CURVE;
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

