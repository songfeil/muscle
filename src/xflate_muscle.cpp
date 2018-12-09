#include "xflate_muscle.h"
#include <eltopo/eltopo.h>
#include <igl/per_vertex_normals.h>
#include <cmath>

 void xflate_muscle(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const int submesh_start,
                    const int submesh_end,
                    const int mode, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew) {
     // Compute desired locations (V + unit_normals!)
    Eigen::MatrixXd N, V_new;
    igl::per_vertex_normals(V, F, N);
     // Only inflate submesh
    Vnew.resizeLike(V);
    for (int i = 0; i < V.rows(); i++){
      Vnew.row(i) = V.row(i);
      if (i >= submesh_start) {
        Vnew.row(i) += N.row(i) * mode * 0.1;
      }
    }
 }

 void xflate_verts(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const std::vector<int> vids,
                    const int mode, // -1 = deflate, 1 = inflate
                    Eigen::MatrixXd & Vnew) {
  Eigen::MatrixXd N, V_new;
  igl::per_vertex_normals(V, F, N);
    // Only inflate submesh
  Vnew = V.replicate(1, 1);
  for (int i : vids) {
    Vnew.row(i) += N.row(i) * mode * 0.01;
  }
}

void xflate_verts_in_xrange(const Eigen::MatrixXd & V,
                    const Eigen::MatrixXi & F,
                    const double xmid,
                    const double xrange,
                    const int mode, // -1 = deflate, 1 = inflate
                    std::set<int> fixed_verts,
                    Eigen::MatrixXd & Vnew){
  Eigen::MatrixXd N, V_new;
  igl::per_vertex_normals(V, F, igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, N);
  Vnew = V.replicate(1, 1);
  for (int i = 0; i < Vnew.rows(); i++) {
    if (fixed_verts.find(i) != fixed_verts.end()) {
      continue;
    }
    double x = Vnew(i, 0);
    double dist = std::abs(x - xmid);
    if (dist < (xrange * 1.5)) {
      //double coeff = - (1.0 / (xrange * xrange)) * (dist * dist) + 1.0;
      double coeff = std::exp((-1.0) * std::pow((x - xmid)/(xrange / 2.0), 4));
      std::cout << coeff << std::endl;
      Vnew.row(i) += N.row(i) * mode * 0.1 * coeff;
    }
  }

  igl::ARAPData data;
  Eigen::VectorXi fixed(s.fixed_vids.size());
  Eigen::MatrixXd V_smooth_fixed(s.fixed_vids.size(), 3);
  int i = 0;
  for (int ind : s.fixed_vids) {
    V_smooth_fixed.row(i) = s.Vm.row(ind);
    fixed(i) = ind;
    i++;
  }
  std::cout << fixed << std::endl;
  igl::arap_precomputation(V_smooth, s.Fm, 3, fixed, data);
  igl::arap_solve(V_smooth_fixed, data, s.Vm);

}


 // EL TOPO BS

 //      // Convert V_new from Eigen matrix to array
//   double *V_newa = new double[3*V_new.rows()];
//   for (int k=0; k<V_new.rows(); k++)
//   {
//     V_newa[3*k] = V_new(k,0);
//     V_newa[3*k+1] = V_new(k,1);
//     V_newa[3*k+2] = V_new(k,2);
//   } 
    
// // Convert V0 from Eigen matrix to array
//   double *Va = new double[3*V.rows()];
//   for (int k=0; k<V.rows(); k++)
//   {
//     Va[3*k] = V(k,0);
//     Va[3*k+1] = V(k,1);
//     Va[3*k+2] = V(k,2);
//   } 
//    std::cout << "Vnew[0]: " << V_newa[0] << std::endl;
//   std::cout << "V[0]: " << Va[0] << std::endl;
//    // Convert F_all from Eigen matrix to array
//   int *Fa = new int[3*F.rows()];
//   for (int k=0; k<F.rows(); k++)
//   {
//     Fa[3*k] = F(k,0);
//     Fa[3*k+1] = F(k,1);
//     Fa[3*k+2] = F(k,2);
//   }
//    // Masses = 1 (all solid)
//   double *masses = new double[V.rows()];
//   for (int i=0; i<V.rows(); i++){
//       masses[i] = 1.0;
//   }
//  // encapsulate all data into an ElTopoMesh
//   ElTopoMesh eltopo_time0;
//   eltopo_time0.num_vertices = V.rows();
//   eltopo_time0.vertex_locations = Va;
//   eltopo_time0.num_triangles = F.rows();
//   eltopo_time0.triangles = Fa;
//   eltopo_time0.vertex_masses = masses;
//      // Set general parameters
//   ElTopoGeneralOptions sim_general_options;
//   // do not print stuff to the console
//   sim_general_options.m_verbose = 0;
//   // do avoid self-intersections
//   sim_general_options.m_collision_safety = 0;
// //    // separation between colliding meshes 
//     //sim_general_options.m_proximity_epsilon = 1e-6; 
//    // Set Simulation parameters
//   ElTopoIntegrationOptions sim_integration_options;
//   sim_integration_options.m_friction_coefficient = 0.0;
//   sim_integration_options.m_dt = 1.0;
//    //
//    double* V_final;
//   double out_dt = 0.0;
//   // We start with 1.0 to step
//   double rest_dt = 1.0;
//   // While we haven't reached final positions
//   ElTopoStaticOperationsOptions sim_static_options;
//   std::cout << "el topooo" << std::endl;
//   while (rest_dt>1e-6)
//   {
//     // ElTopoMesh new_mesh;
//     // ElTopoDefragInformation defrag;
//     // sim_static_options.m_allow_non_manifold = false;
//     // sim_static_options.m_allow_topology_changes = false;
//     // call Eltopo main function
//     el_topo_integrate(&eltopo_time0, V_newa, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
//     //el_topo_static_operations(&eltopo_time0, &sim_general_options, &sim_static_options, &defrag, &eltopo_time0);
//     std::cout << "out_dt = " << out_dt << std::endl;
//     // update the rest to go
//     rest_dt = (1-out_dt)*rest_dt;
//     // if we haven't reached final positions, print how much we have stepped
//     // and update vertex positions 
//     if (out_dt < 1.0)
//     {
//       std::cout << "out_dt = " << out_dt << std::endl;
//       eltopo_time0.vertex_locations = V_final;
//     }
//   }
//    std::cout << "V_final[0]: " << V_final[0] << std::endl;
//   // output corrected velocities
//   Vnew.resize(submesh_end - submesh_start, 3);
//   for (int k=submesh_start; k<submesh_end; k++)
//   {
//     Vnew(k - submesh_start, 0) = V_final[3*k];
//     Vnew(k - submesh_start, 1) = V_final[3*k+1];
//     Vnew(k - submesh_start, 2) = V_final[3*k+2];
//   }