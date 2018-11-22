#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_ray.h>
#include <igl/combine.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/read_triangle_mesh.h>

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <vector>
#include <iostream>

#include "generate_bone.h"
#include "generate_muscle.h"
#include "ray_intersect_plane.h"
//#include "inflate_muscle.h"

/* Input mode enum */
enum Mode
{
  NONE = 0,
  BONE = 1,
  MUSCLE = 2,
  FACE_SELECT = 3,
};

/* Color consts */
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);

/* GLOBAL STATE :-) */

// Mesh
std::vector<Eigen::MatrixXd> VV;
std::vector<Eigen::MatrixXi> FF;
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd face_colors;
std::set<int> selected_faces;

// Interface
Mode mode = NONE;
Eigen::Vector3d last_mouse;
Eigen::MatrixXd control_points;
int n_points = 0;

Eigen::MatrixXd Vbone;
Eigen::MatrixXi Fbone;
Eigen::Vector3d bp0;


int main(int argc, char *argv[])
{

  /* Instantiate libigl Viewer */
  igl::opengl::glfw::Viewer viewer;
  std::cout << "========= MUSCLE/BONE GENERATION =========" << std::endl;
  std::cout << " B,b \t Bone input mode \n M,m \t Muscle input mode \n[space]  No mode (can drag around scene again)" << std::endl;

  // Import our template bone mesh
  igl::read_triangle_mesh("../data/bone.obj",Vbone,Fbone);
  
  /* Define update functions for viewer */
  const auto & update = [&]()
  {

    /* If FACE_SELECT, only update colours */
    if (mode == FACE_SELECT) {
      viewer.data().set_colors(face_colors);
    }

    /* If BONE or MUSCLE generated, create a new combination of meshes */
    // Setting multiple meshes
    // http://www.alecjacobson.com/weblog/?p=4679
    else {
      
      viewer.data().clear();
      viewer.data().set_points(control_points, mode == BONE? orange : blue);

      if (VV.size() > 0) {

        // New mesh!
        igl::combine(VV,FF,V,F);
        viewer.data().set_mesh(V,F);

        // Reset colors
        face_colors.resize(F.rows(), 3);
        for (int i = 0; i < F.rows(); i++) {
          face_colors.row(i) = green;
        }

        viewer.data().set_colors(face_colors);
      }
    }
  };

  /* Mouse down */
  viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {

    /* Not NONE */
    if (mode > 0) {
      last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);

      /* FACE_SELECT */
      if (mode == FACE_SELECT) {
        // Ray cast to find closest face, color it.
        int fid;
        Eigen::Vector3f bary;
        Eigen::RowVector2f mouse;
        mouse << last_mouse.head(2)(0), last_mouse.head(2)(1);

        if (igl::unproject_onto_mesh(mouse, viewer.core.view, viewer.core.proj, viewer.core.viewport, V, F, fid, bary)) {
          selected_faces.insert(fid);
          face_colors.row(fid) = orange;
        }
      }

      /* Point input modes */
      else {
        Eigen::Vector3d source, dir, intersection;
        Eigen::Vector2d p;
        p << last_mouse.head(2)(0), last_mouse.head(2)(1);
        igl::unproject_ray(p, viewer.core.view, viewer.core.proj, viewer.core.viewport, source, dir);
        ray_intersect_plane(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), source, dir, intersection);

        /* BONE */
        if (mode == BONE) {
          
          if (n_points == 0) {
            Eigen::MatrixXd Vnew = Vbone.replicate(1, 1);
            start_flinstone_bone(intersection, Vnew);
            VV.push_back(Vnew);
            FF.push_back(Fbone);
            bp0 = intersection;
            n_points++;
          } 
          else if (n_points == 1) {
            Eigen::MatrixXd Vnew = Vbone.replicate(1, 1);
            transform_flinstone_bone(bp0, intersection, Vnew);
            VV.pop_back();
            VV.push_back(Vnew);
            n_points = 0;
          }
          control_points.conservativeResize(control_points.rows() + 1, 3);
          control_points.row(control_points.rows() - 1) = intersection;

        }

        /* MUSCLE */
        else if (mode == MUSCLE) {

          control_points.conservativeResize(control_points.rows() + 1, 3);
          control_points.row(control_points.rows() - 1) = intersection;
          n_points++;

          // For every 3 points, generate a muscle!
          if (n_points == 3) {
            generate_muscle(control_points, n_points, V, F, selected_faces, VV, FF);
            n_points = 0;
          }
        }
      }

        update();
        return true;

    }

    return false;

  };

  viewer.callback_mouse_move = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    if ((mode == BONE) && (n_points == 1)) {
      last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);
      Eigen::Vector3d source, dir, intersection;
      Eigen::Vector2d p;
      p << last_mouse.head(2)(0), last_mouse.head(2)(1);
      igl::unproject_ray(p, viewer.core.view, viewer.core.proj, viewer.core.viewport, source, dir);
      ray_intersect_plane(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), source, dir, intersection);

      Eigen::MatrixXd Vnew = Vbone.replicate(1, 1);
      transform_flinstone_bone(bp0, intersection, Vnew);
      VV.pop_back();
      VV.push_back(Vnew);
      update();
      return true;
    }
    return false;
  };

  /* Key presses */
  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch(key) {
      case 'B':
      case 'b':
      {
         mode = BONE;
         control_points.resize(0, 3);
         n_points = 0;
        break;
      }
      case 'M':
      case 'm':
      {
         mode = MUSCLE;
         control_points.resize(0, 3);
         n_points = 0;
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
          generate_bones(control_points, VV, FF);
        }
        break;
      }
      case 'i':
      {
        Eigen::MatrixXd Vnew;
        Eigen::MatrixXi Fnew;
        Eigen::MatrixXd muscle_mesh = VV.back(); // Assuming muscle was last mesh added
        int start = V.rows() - muscle_mesh.rows();
        int end = V.rows();
        //inflate_muscle(V, F, start, end, Vnew);
        VV.pop_back();
        VV.push_back(Vnew);
        std::cout << "INFLATE" << std::endl;
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

