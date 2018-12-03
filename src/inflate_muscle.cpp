// #include "inflate_muscle.h"
// #include <eltopo/eltopo.h>
// #include <eltopo/subdivisionscheme.h>
// #include <igl/per_vertex_normals.h>

// Eigen::VectorXd dual_volumes(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F) {

//     Eigen::VectorXd Vols = Eigen::VectorXd::Zero(V.rows());

//     auto double_vol = [&](int ai, int bi, int ci) -> double {
//         Eigen::Vector3d a = V.row(ai);
//         Eigen::Vector3d b = V.row(bi);
//         Eigen::Vector3d c = V.row(ci);
//         return (b-a).cross(c-a).norm();
//     };
//     for(int i = 0; i < F.rows(); ++i) {
//         auto f = F.row(i);
//         double v = double_vol(f(0),f(1),f(2));
//         for(int j = 0; j < 3; ++j) {
//             Vols(f(j)) += v;
//         }
//     }
//     Vols /= 6;// a third for each cell, half for the double-vols
//     return Vols;
// }

// void inflate_muscle(const Eigen::MatrixXd & V,
//                    const Eigen::MatrixXi & F,
//                    int submesh_start,
//                    int submesh_end,
//                    Eigen::MatrixXd & Vnew) {

//    // Compute desired locations (V + unit_normals!)
//    Eigen::MatrixXd N, V_new;
//    igl::per_vertex_normals(V, F, N);

//    // Only inflate submesh
//    V_new.resizeLike(V);
//    for (int i = 0; i < V.rows(); i++){
//      V_new.row(i) = V.row(i);
//      if (i >= submesh_start) {
//        V_new.row(i) += N.row(i) * 0.1;
//      }
//    }

//     //ElTopoTracker tracker(V, F);
//     // for (int i = 0; i < 10; i++){
//     //     tracker.step(V_new, 0.1);
//     // }

//    // Convert V_new from Eigen matrix to array
//  double *V_newa = new double[3*V_new.rows()];
//  for (int k=0; k<V_new.rows(); k++)
//  {
//    V_newa[3*k] = V_new(k,0);
//    V_newa[3*k+1] = V_new(k,1);
//    V_newa[3*k+2] = V_new(k,2);
//  }

// // Convert V0 from Eigen matrix to array
//  double *Va = new double[3*V.rows()];
//  for (int k=0; k<V.rows(); k++)
//  {
//    Va[3*k] = V(k,0);
//    Va[3*k+1] = V(k,1);
//    Va[3*k+2] = V(k,2);
//  }

//  std::cout << "Vnew[0]: " << V_newa[0] << std::endl;
//  std::cout << "V[0]: " << Va[0] << std::endl;

//  // Convert F_all from Eigen matrix to array
//  int *Fa = new int[3*F.rows()];
//  for (int k=0; k<F.rows(); k++)
//  {
//    Fa[3*k] = F(k,0);
//    Fa[3*k+1] = F(k,1);
//    Fa[3*k+2] = F(k,2);
//  }

//  // thanks mike?
//  Eigen::VectorXd vols = dual_volumes(V, F);
//  double *masses = new double[V.rows()];
//  for (int i=0; i<vols.rows(); i++){
//      masses[i] = 1;
//  }

// // encapsulate all data into an ElTopoMesh
//  ElTopoMesh eltopo_time0;
//  eltopo_time0.num_vertices = V.rows();
//  eltopo_time0.vertex_locations = Va;
//  eltopo_time0.num_triangles = F.rows();
//  eltopo_time0.triangles = Fa;
//  eltopo_time0.vertex_masses = masses;

//    // Set general parameters
//  ElTopoGeneralOptions sim_general_options;
//  // do not print stuff to the console

//  // Set Simulation parameters
//  ElTopoIntegrationOptions sim_integration_options;

//   ElTopoStaticOperationsOptions sim_static_op_options;
//   sim_static_op_options.m_subdivision_scheme = new ButterflyScheme();
//   sim_static_op_options.m_allow_non_manifold = false;

//   ElTopoDefragInformation defrag_info;

//  //
//  double* V_final;
//  double out_dt = 0.0;
//  // We start with 1.0 to step
//  double rest_dt = 1.0;
//  // While we haven't reached final positions
//  std::cout << "el topooo" << std::endl;
//  Vnew.resize(submesh_end - submesh_start, 3);
//  while (rest_dt>1e-6)
//  {
//    // call Eltopo main function
//    //el_topo_static_operations(&eltopo_time0, &sim_general_options, &sim_static_op_options, &defrag_info, &eltopo_time0);
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

//  std::cout << "V_final[0]: " << V_final[0] << std::endl;
//  // output corrected velocities
//  Vnew.resize(submesh_end - submesh_start, 3);
//  for (int k=submesh_start; k<submesh_end; k++)
//  {
//    Vnew(k - submesh_start, 0) = V_final[3*k];
//    Vnew(k - submesh_start, 1) = V_final[3*k+1];
//    Vnew(k - submesh_start, 2) = V_final[3*k+2];
//  }

// }