//#include "inflate_muscle.h"
//#include <eltopo/eltopo.h>
//#include <igl/per_vertex_normals.h>
//
//void inflate_muscle(const Eigen::MatrixXd & V,
//                    const Eigen::MatrixXi & F,
//                    int submesh_start,
//                    int submesh_end,
//                    Eigen::MatrixXd & Vnew) {
//
//    // Compute desired locations (V + unit_normals!)
//    Eigen::MatrixXd N, V_new;
//    igl::per_vertex_normals(V, F, N);
//
//    // Only inflate submesh
//    V_new.resizeLike(V);
//    for (int i = 0; i < V.rows(); i++){
//      V_new.row(i) = V.row(i);
//      if (i >= submesh_start) {
//        V_new.row(i) += N.row(i) * 0.1;
//      }
//    }
//
//    // Convert V_new from Eigen matrix to array
//  double *V_newa = new double[3*V_new.rows()];
//  for (int k=0; k<V_new.rows(); k++)
//  {
//    V_newa[3*k] = V_new(k,0);
//    V_newa[3*k+1] = V_new(k,1);
//    V_newa[3*k+2] = V_new(k,2);
//  }
//
//// Convert V0 from Eigen matrix to array
//  double *Va = new double[3*V.rows()];
//  for (int k=0; k<V.rows(); k++)
//  {
//    Va[3*k] = V(k,0);
//    Va[3*k+1] = V(k,1);
//    Va[3*k+2] = V(k,2);
//  }
//
//  std::cout << "Vnew[0]: " << V_newa[0] << std::endl;
//  std::cout << "V[0]: " << Va[0] << std::endl;
//
//  // Convert F_all from Eigen matrix to array
//  int *Fa = new int[3*F.rows()];
//  for (int k=0; k<F.rows(); k++)
//  {
//    Fa[3*k] = F(k,0);
//    Fa[3*k+1] = F(k,1);
//    Fa[3*k+2] = F(k,2);
//  }
//
//  // Masses = 1 (all solid)
//  double *masses = new double[V.rows()];
//  for (int i=0; i<V.rows(); i++){
//      masses[i] = 1.0;
//  }
//
//// encapsulate all data into an ElTopoMesh
//  ElTopoMesh eltopo_time0;
//  eltopo_time0.num_vertices = V.rows();
//  eltopo_time0.vertex_locations = Va;
//  eltopo_time0.num_triangles = F.rows();
//  eltopo_time0.triangles = Fa;
//  eltopo_time0.vertex_masses = masses;
//
//    // Set general parameters
//  ElTopoGeneralOptions sim_general_options;
//  // do not print stuff to the console
//  sim_general_options.m_verbose = 0;
//  // do avoid self-intersections
//  sim_general_options.m_collision_safety = 0;
////    // separation between colliding meshes
//    //sim_general_options.m_proximity_epsilon = 1e-6;
//
//
//  // Set Simulation parameters
//  ElTopoIntegrationOptions sim_integration_options;
//  sim_integration_options.m_friction_coefficient = 0.0;
//  sim_integration_options.m_dt = 1.0;
//
//  //
//
//
//  double* V_final;
//  double out_dt = 0.0;
//  // We start with 1.0 to step
//  double rest_dt = 1.0;
//  // While we haven't reached final positions
//  std::cout << "el topooo" << std::endl;
//  while (rest_dt>1e-6)
//  {
//    // call Eltopo main function
//    el_topo_integrate(&eltopo_time0, V_newa, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
//    std::cout << "out_dt = " << out_dt << std::endl;
//    // update the rest to go
//    rest_dt = (1-out_dt)*rest_dt;
//    // if we haven't reached final positions, print how much we have stepped
//    // and update vertex positions
//    if (out_dt < 1.0)
//    {
//      std::cout << "out_dt = " << out_dt << std::endl;
//      eltopo_time0.vertex_locations = V_final;
//    }
//  }
//
//  std::cout << "V_final[0]: " << V_final[0] << std::endl;
//  // output corrected velocities
//  Vnew.resize(submesh_end - submesh_start, 3);
//  for (int k=submesh_start; k<submesh_end; k++)
//  {
//    Vnew(k - submesh_start, 0) = V_final[3*k];
//    Vnew(k - submesh_start, 1) = V_final[3*k+1];
//    Vnew(k - submesh_start, 2) = V_final[3*k+2];
//  }
//
//}