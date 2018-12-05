#include <igl/opengl/glfw/Viewer.h>
#include <igl/combine.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <igl/adjacency_matrix.h>

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <vector>
#include <string>
#include <iostream>
#include <stack>

#include "ui_helpers.h"
#include "generate_bone.h"
#include "generate_muscle.h"
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
const Eigen::RowVector3d red(1.0,0.0,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);
Eigen::RowVector3d patch_colors[3] = {yellow, red, blue};

/* GLOBAL STATE :-) */

// Undoable
struct State
{
  int id = 0;
  // Lists of meshes
  std::vector<Eigen::MatrixXd> VV; 
  std::vector<Eigen::MatrixXi> FF;
  // Total combination of all meshes to be displayed
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Selected faces
  Eigen::MatrixXd face_colors;
  std::vector<std::vector<int>> patch_fids; // for coloring
  std::vector<Eigen::MatrixXi> patch_faces;
  // Muscle control points
  Eigen::MatrixXd control_points;
  int n_points = 0;
  // Bone start point
  Eigen::Vector3d bp0;
} s;

// Non-undoable state
Mode mode = NONE;
Eigen::Vector3d last_mouse;
bool bone_started = false;
bool mouse_down = false;
int hover_point_index = -1;

// Bone mesh template
Eigen::MatrixXd Vbone;
Eigen::MatrixXi Fbone;

int main(int argc, char *argv[])
{

  // Undo Management
  std::stack<State> undo_stack,redo_stack;
  const auto push_undo = [&](State & _s=s)
  {
    undo_stack.push(_s);
    // clear
    redo_stack = std::stack<State>();
  };
  const auto undo = [&]()
  {
    if(!undo_stack.empty())
    {
      redo_stack.push(s);
      s = undo_stack.top();
      undo_stack.pop();
    }
  };
  const auto redo = [&]()
  {
    if(!redo_stack.empty())
    {
      undo_stack.push(s);
      s = redo_stack.top();
      redo_stack.pop();
    }
  };

  /* Instantiate libigl Viewer */
  igl::opengl::glfw::Viewer viewer;
  std::cout << "========= MUSCLE/BONE GENERATION =========" << std::endl;
  std::cout << " B,b \t Bone input mode \n M,m \t Muscle input mode \n[space]  No mode (can drag around scene again)" << std::endl;

  // Import our template bone mesh
  igl::read_triangle_mesh("../data/bone.obj",Vbone,Fbone);
  
  /* Define update functions for viewer */
  const auto & update = [&](bool update_mesh)
  {
      
    viewer.data().clear();
    viewer.data().set_points(s.control_points, yellow);

    if (s.VV.size() > 0) {
      
      // New mesh!
      if (update_mesh) {
        igl::combine(s.VV,s.FF,s.V,s.F);
        viewer.data().set_mesh(s.V,s.F);
      }
      
      // Reset colors
      s.face_colors.resize(s.F.rows(), 3);
      for (int i = 0; i < s.F.rows(); i++) {
        s.face_colors.row(i) = green;
      }

      // Color each patch differently
      int patch_num = 0;
      for (const auto & patch : s.patch_fids) {
        int color = patch_num % 3;
        for (int fid : patch) {
          s.face_colors.row(fid) = patch_colors[color];
        }
        patch_num++;
      }

      viewer.data().set_colors(s.face_colors);
    }

  };

  /* Mouse down */
  viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    mouse_down = true;
    last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);

    /* Not NONE */
    if (mode > 0) {

      /* FACE_SELECT */
      if (mode == FACE_SELECT) {
        // Ray cast to find closest face, color it.
        int fid;
        Eigen::Vector3f bary;
        Eigen::RowVector2f mouse;
        mouse << last_mouse.head(2)(0), last_mouse.head(2)(1);

        if (igl::unproject_onto_mesh(mouse, viewer.core.view, viewer.core.proj, viewer.core.viewport, s.V, s.F, fid, bary)) {

          push_undo(s);
          add_face_to_patch(fid, s.F.row(fid), s.patch_faces, s.patch_fids);
          update(true);
          return true;
          
        }

      }

      /* Point input modes */
      else {
        Eigen::Vector3d intersection;
        intersection_with_xy_plane(viewer, last_mouse, intersection);

        /* BONE */
        if (mode == BONE) {
          
          // Start of bone
          if (!bone_started) {

            push_undo(s);
            
            Eigen::MatrixXd Vnew = Vbone.replicate(1, 1);
            start_flinstone_bone(intersection, Vnew);
            s.VV.push_back(Vnew);
            s.FF.push_back(Fbone);
            s.bp0 = intersection;
            bone_started = true;
          } 
          // End of bone
          else {

            Eigen::MatrixXd Vnew = Vbone.replicate(1, 1);
            transform_flinstone_bone(s.bp0, intersection, Vnew);
            s.VV.pop_back();
            s.VV.push_back(Vnew);
            bone_started = false;

          }

          update(true);

        }

        /* MUSCLE */
        else if (mode == MUSCLE) {

          push_undo(s);

          s.control_points.conservativeResize(s.control_points.rows() + 1, 3);
          s.control_points.row(s.control_points.rows() - 1) = intersection;
          s.n_points++;

          update(true);

        }
      }
      return true;

    }

    return false;

  };

  viewer.callback_mouse_move = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
  {
    last_mouse = Eigen::Vector3d(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);
    
    // Wiggle around bone after one point has been put down!
    if ((mode == BONE) && (bone_started)) {
      
      Eigen::Vector3d intersection;
      intersection_with_xy_plane(viewer, last_mouse, intersection);

      Eigen::MatrixXd Vnew = Vbone.replicate(1, 1);
      transform_flinstone_bone(s.bp0, intersection, Vnew);
      s.VV.pop_back();
      s.VV.push_back(Vnew);

      update(true);
      return true;
    }
  
    else if (mode == MUSCLE) {

      Eigen::Vector3d intersection;
      intersection_with_xy_plane(viewer, last_mouse, intersection);

      // Continuous "sketching" input for muscle points
      if (mouse_down) {

        push_undo(s);

        s.control_points.conservativeResize(s.control_points.rows() + 1, 3);
        s.control_points.row(s.control_points.rows() - 1) = intersection;
        s.n_points++;

        update(true);
        return true;
      }

      // // Hover existing point + click to delete
      // else {

      //   for (int i = 0; i < s.control_points.rows(); i++) {
      //     auto p = s.control_points.row(i);
      //     if ((intersection - p).norm() < viewer.data().point_size) {

      //     }
      //   }

      // }
      
    }
    return false;
  };

  /* Mouse up */
  viewer.callback_mouse_up = 
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
    {
      mouse_down = false;
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
        // Bone input mode
        mode = BONE;
        break;
      }
      case 'M':
      case 'm':
      {
        // Muscle input mode
        if (mode == MUSCLE) { // click M again to clear? lol
          s.control_points.resize(0, 3);
          s.n_points = 0;
          update(true);
        } else {
          mode = MUSCLE;
        }
        break;
      }
      case ' ':
      {
        // Panning/scene browsing
        mode = NONE;
        break;
      }
      case 'f':
      {
        // Face selecting mode
        mode = FACE_SELECT;
        break;
      }
      case 'G':
      case 'g':
      {
        // Generate muscle from points
        if (s.n_points > 2) { // Need some points to work with...
          // NO LONGER HAVE S.SELECTED_FACES...
          // Just picking a face from each patch......
          // Temporary!
          std::set<int> face_per_patch; 
          for (const auto & patch : s.patch_fids) {
            face_per_patch.insert(patch.at(0));
          }
//          generate_muscle(s.control_points, s.n_points, s.V, s.F, s.patch_faces, s.VV, s.FF);
          generate_muscle_multiface(s.control_points, s.n_points, s.V, s.F, s.patch_faces, s.VV, s.FF);
          s.control_points.resize(0, 3);
          s.n_points = 0;

          update(true);
        }
        break;
      }
      case 'i':
      {
        // // Inflate muscle mesh
        // Eigen::MatrixXd Vnew;
        // Eigen::MatrixXi Fnew;
        // Eigen::MatrixXd muscle_mesh = s.VV.back(); // Assuming muscle was last mesh added
        // int start = s.V.rows() - muscle_mesh.rows();
        // int end = s.V.rows();
        // inflate_muscle(s.V, s.F, start, end, Vnew);
        // s.VV.pop_back();
        // s.VV.push_back(Vnew);
        // std::cout << "INFLATE" << std::endl;
        // break;
      }
      case 's':
      {
        Eigen::MatrixXd empty;
        igl::writeOBJ("output_all.obj", s.V, s.F, empty, empty, empty, empty);
        for (int i = 0; i < s.VV.size(); i++) {
          std::string suffix = ".obj";
          std::string filename = std::to_string(i) + suffix;
          igl::writeOBJ(filename, s.VV[i], s.FF[i], empty, empty, empty, empty);
        }
        break;
      }
    }
    return true;
  };

  // Special callback for handling undo
  viewer.callback_key_down = 
    [&](igl::opengl::glfw::Viewer &, unsigned char key, int mod)->bool
  {
    if(key == 'Z' && (mod & GLFW_MOD_SUPER))
    {
      std::cout << "undo/redo?" << std::endl;
      (mod & GLFW_MOD_SHIFT) ? redo() : undo();
      update(true);
      return true;
    }
    return false;
  };

  // Launch a viewer instance
  viewer.launch();
  return 0;
}

