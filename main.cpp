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
#include <iostream>

int main(int argc, char *argv[])
{
  std::vector<Eigen::MatrixXd> VV;
  std::vector<Eigen::MatrixXi> FF;
  // Load in a mesh
//  igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/knight.off", V, F);

//
//  Eigen::MatrixXi E = edges(F);
//  int Chi = euler_characteristic(F);
//  std::cout<<"Edge list E is "<<E.rows()<<"x"<<E.cols()<<std::endl;
//  std::cout<<"Euler Characteristic: "<<Chi<<std::endl;

  // OUR CODE HERE!
  //generate_bone(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 3, 5), V, F);

  Eigen::MatrixXd B, N;
  bezier(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 2.5), Eigen::Vector3d(0, 0, 5), 10, B, N);

//  std::cout << B << std::endl;
  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  cylinder(12, 22, CV, CF);
  VV.push_back(CV);
  FF.push_back(CF);

  // Create a libigl Viewer object 
  igl::opengl::glfw::Viewer viewer;

  // Color consts
  const Eigen::RowVector3d orange(1.0,0.7,0.2);
  const Eigen::RowVector3d yellow(1.0,0.9,0.2);
  const Eigen::RowVector3d blue(0.2,0.3,0.8);
  const Eigen::RowVector3d green(0.2,0.6,0.3);

  // User input data
  Eigen::Vector3d last_mouse;
  Eigen::MatrixXd bone_points;
  int n_bp = 0;
  bool creating_bones = true; //default

  // Define update functions for viewer
  const auto & update = [&]()
  {
    if (creating_bones) {
      viewer.data().clear();
      viewer.data().set_points(bone_points,orange);
      // Setting multiple meshes
      // http://www.alecjacobson.com/weblog/?p=4679
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
    last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);
    if (creating_bones && (n_bp < 2)) {
      // Ray cast to x-y plane
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
      // Add intersection to list of bone points
      bone_points.conservativeResize(bone_points.rows() + 1, 3);
      bone_points.row(bone_points.rows() - 1) = intersection;
      n_bp++;
      // For every 2 points, generate a bone
      if (n_bp == 2){
//        std::cout << "GEN BONE!" << std::endl;
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::Vector3d p1 = bone_points.row(bone_points.rows() - 1);
        Eigen::Vector3d p2 = bone_points.row(bone_points.rows() - 2);
        generate_bone(p1, p2, V, F);
//        std::cout << "all bone points: \n" << bone_points << std::endl;
//        std::cout << "new V: \n" << V << std::endl;
        // Add to list of meshes
        VV.push_back(V);
        FF.push_back(F);
//        std::cout << "num meshes: \n" << VV.size() << std::endl;
        n_bp = 0;
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
        creating_bones = true;
        break;
      }
      case ' ':
      {
        creating_bones = false;
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

