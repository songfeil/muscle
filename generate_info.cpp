#include <igl/AABB.h>
#include <igl/signed_distance.h>
#include <igl/readOBJ.h>
#include <igl/readMSH.h>
#include <vector>
#include <iostream>

typedef Eigen::MatrixXd VMAT;
typedef Eigen::MatrixXi TMAT;
typedef Eigen::MatrixXi FMAT;
typedef Eigen::RowVector4i Tet;
typedef Eigen::RowVector3d Vertex;

class info {
 public:
  std::vector<int> bone1Idx;
  std::vector<int> bone2Idx;
  std::vector<int> muscleIdx;
};

info generate_info(
      const VMAT & Vb1,
      const TMAT & Fb1,
      const VMAT & Vb2,
      const TMAT & Fb2,
      const VMAT & Vm,
      const TMAT & Fm,
      const VMAT & V_all,
      const TMAT & T_all
    ) {

  info result;

  Eigen::MatrixXd Q;
  Eigen::VectorXd Sb1, Sb2, Sm, pH0;
  Eigen::MatrixXd pH1, pH2;

  // Calculate central point
  Q.resize(T_all.rows(), 3);
  for (int i = 0; i < T_all.rows(); i++) {
    Tet t = T_all.row(i);
    Q.row(i) = (V_all.row(t(0)) + V_all.row(t(1)) + V_all.row(t(2)) + V_all.row(t(3))) * 0.25;
  }

  // Signed distance
  igl::SignedDistanceType type = igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
  igl::signed_distance(Q, Vb1, Fb1, type, Sb1, pH0, pH1, pH2);
  igl::signed_distance(Q, Vb2, Fb2, type, Sb2, pH0, pH1, pH2);
  igl::signed_distance(Q, Vm, Fm, type, Sm, pH0, pH1, pH2);


  for (int i = 0; i < T_all.rows(); i++) {
    if (Sb1(i) < 0) {
      // Found in bone1
      result.bone1Idx.push_back(i);
    } else if (Sb2(i) < 0) {
      result.bone2Idx.push_back(i);
    } else if (Sm(i) < 0) {
      result.muscleIdx.push_back(i);
    } else {
      std::cout << "Found point #" << i << " not found in any mesh." << std::endl;
      double db1, db2, dm;
      db1 = Sb1(i);
      db2 = Sb2(i);
      dm = Sm(i);

      if (db1 <= db2 && db1 <= dm) {
        result.bone1Idx.push_back(i);
      } else if (db2 <= db1 && db2 <= dm) {
        result.bone2Idx.push_back(i);
      } else {
        result.muscleIdx.push_back(i);
      }
    }
  }

  return result;
};

int main(int argc, char ** argv) {
  VMAT Vb1, Vb2, Vm, V_all;
  FMAT Fb1, Fb2, Fm;
  TMAT T_all;
  igl::readOBJ(argv[1], Vb1, Fb1);
  igl::readOBJ(argv[2], Vb2, Fb2);
  igl::readOBJ(argv[3], Vm, Fm);
  igl::readMSH(argv[4], V_all, T_all);

  info result = generate_info(Vb1, Fb1, Vb2, Fb2, Vm, Fm, V_all, T_all);
  std::cout << "Bone 1 has " << result.bone1Idx.size() << " tets" << std::endl;
  std::cout << "Bone 2 has " << result.bone2Idx.size() << " tets" << std::endl;
  std::cout << "Muscle has " << result.muscleIdx.size() << " tets" << std::endl;

  return 0;
}