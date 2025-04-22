#include "permutation.h"
#include <iostream>

void rotY(VectorXf& x, int i, float angle) {
  const float cs = std::cos(2.0f * float(EIGEN_PI) * angle / 360.0f);
  const float sn = std::sin(2.0f * float(EIGEN_PI) * angle / 360.0f);
  Matrix3f m;
  m << cs, 0.0f, -sn,
    0.0f, 1.0f, 0.0f,
    sn, 0.0f, cs;
  Eigen::Map<Vector3f>(x.data() + i * 3) = (m * Eigen::Map<Vector3f>(x.data() + i * 3)).eval();
}
void rotZ(VectorXf& x, int i, float angle) {
  const float cs = std::cos(2.0f * float(EIGEN_PI) * angle / 360.0f);
  const float sn = std::sin(2.0f * float(EIGEN_PI) * angle / 360.0f);
  Matrix3f m;
  m << cs, -sn, 0.0f,
    sn, cs, 0.0f,
    0.0f, 0.0f, 1.0f;
  Eigen::Map<Vector3f>(x.data() + i * 3) = (m * Eigen::Map<Vector3f>(x.data() + i * 3)).eval();
}

//0 DOF
void set_000(VectorXf& x, int i) {
  x[i * 3 + 0] = 0.0f;
  x[i * 3 + 1] = 0.0f;
  x[i * 3 + 2] = 0.0f;
}
//1 DOF
void set_x00(VectorXf& x, int i) {
  x[i * 3 + 1] = 0.0f;
  x[i * 3 + 2] = 0.0f;
}
void set_0y0(VectorXf& x, int i) {
  x[i * 3 + 0] = 0.0f;
  x[i * 3 + 2] = 0.0f;
}
void set_00z(VectorXf& x, int i) {
  x[i * 3 + 0] = 0.0f;
  x[i * 3 + 1] = 0.0f;
}
void set_xx0(VectorXf& x, int i) {
  x[i * 3 + 1] = x[i * 3 + 0];
  x[i * 3 + 2] = 0.0f;
}
void set_xn0(VectorXf& x, int i) {
  x[i * 3 + 1] = -x[i * 3 + 0];
  x[i * 3 + 2] = 0.0f;
}
void set_x0x(VectorXf& x, int i) {
  x[i * 3 + 1] = 0.0f;
  x[i * 3 + 2] = x[i * 3 + 0];
}
void set_x0n(VectorXf& x, int i) {
  x[i * 3 + 1] = 0.0f;
  x[i * 3 + 2] = -x[i * 3 + 0];
}
void set_0yy(VectorXf& x, int i) {
  x[i * 3 + 0] = 0.0f;
  x[i * 3 + 2] = x[i * 3 + 1];
}
void set_0yn(VectorXf& x, int i) {
  x[i * 3 + 0] = 0.0f;
  x[i * 3 + 2] = -x[i * 3 + 1];
}
void set_xxx(VectorXf& x, int i) {
  x[i * 3 + 1] = x[i * 3 + 0];
  x[i * 3 + 2] = x[i * 3 + 0];
}
void set_xnx(VectorXf& x, int i) {
  x[i * 3 + 1] = -x[i * 3 + 0];
  x[i * 3 + 2] = x[i * 3 + 0];
}
void set_xxn(VectorXf& x, int i) {
  x[i * 3 + 1] = x[i * 3 + 0];
  x[i * 3 + 2] = -x[i * 3 + 0];
}
void set_xnn(VectorXf& x, int i) {
  x[i * 3 + 1] = -x[i * 3 + 0];
  x[i * 3 + 2] = -x[i * 3 + 0];
}
//2 DOF
void set_0yz(VectorXf& x, int i) {
  x[i * 3 + 0] = 0.0f;
}
void set_x0z(VectorXf& x, int i) {
  x[i * 3 + 1] = 0.0f;
}
void set_xy0(VectorXf& x, int i) {
  x[i * 3 + 2] = 0.0f;
}

void x_symmetric(VectorXf& x, int i, int j) {
  x[j * 3 + 0] = x[i * 3 + 0];
  x[j * 3 + 1] = -x[i * 3 + 1];
  x[j * 3 + 2] = -x[i * 3 + 2];
}
void z_symmetric(VectorXf& x, int i, int j) {
  x[j * 3 + 0] = -x[i * 3 + 0];
  x[j * 3 + 1] = -x[i * 3 + 1];
  x[j * 3 + 2] = x[i * 3 + 2];
}

void pr_symmetric(VectorXf& x, int i, int j) {
  x[j * 3 + 0] = -x[i * 3 + 0];
  x[j * 3 + 1] = -x[i * 3 + 1];
  x[j * 3 + 2] = -x[i * 3 + 2];
}
void xref_symmetric(VectorXf& x, int i, int j) {
  x[j * 3 + 0] = -x[i * 3 + 0];
  x[j * 3 + 1] = x[i * 3 + 1];
  x[j * 3 + 2] = x[i * 3 + 2];
}
void zref_symmetric(VectorXf& x, int i, int j) {
  x[j * 3 + 0] = x[i * 3 + 0];
  x[j * 3 + 1] = x[i * 3 + 1];
  x[j * 3 + 2] = -x[i * 3 + 2];
}

void z_offaxis_sym3(VectorXf& x, int i, int j, int k) {
  x[j * 3 + 0] = x[i * 3 + 1];
  x[j * 3 + 1] = x[i * 3 + 2];
  x[j * 3 + 2] = x[i * 3 + 0];
  x[k * 3 + 0] = x[i * 3 + 2];
  x[k * 3 + 1] = x[i * 3 + 0];
  x[k * 3 + 2] = x[i * 3 + 1];
}

void z_dihedral2(VectorXf& x, const Face& s) {
  x[s[1] * 3 + 0] = -x[s[0] * 3 + 0];
  x[s[1] * 3 + 1] = -x[s[0] * 3 + 1];
  x[s[1] * 3 + 2] = x[s[0] * 3 + 2];
  x[s[2] * 3 + 0] = x[s[0] * 3 + 0];
  x[s[2] * 3 + 1] = -x[s[0] * 3 + 1];
  x[s[2] * 3 + 2] = -x[s[0] * 3 + 2];
  x[s[3] * 3 + 0] = -x[s[0] * 3 + 0];
  x[s[3] * 3 + 1] = x[s[0] * 3 + 1];
  x[s[3] * 3 + 2] = -x[s[0] * 3 + 2];
}

void z_symmetric(VectorXf& x, const Face& s, bool alternate = false, int quotient=1) {
  const int size = (int)s.size();
  const float a = (alternate ? -1.0f : 1.0f);
  const float cs = std::cos(quotient * 2.0f * float(EIGEN_PI) / size);
  const float sn = std::sin(quotient * 2.0f * float(EIGEN_PI) / size);
  Matrix3f m;
  m << cs, -sn, 0.0f,
       sn, cs, 0.0f,
       0.0f, 0.0f, a;
  for (int i = 0; i < size - quotient; i += quotient) {
    const int prev_ix = i;
    const int next_ix = (i + quotient) % size;
    Eigen::Map<Vector3f>(x.data() + s[next_ix] * 3) = m * Eigen::Map<Vector3f>(x.data() + s[prev_ix] * 3);
  }
}
void z_dihedral(VectorXf& x, const Face& s) {
  const int size = (int)s.size();
  const float cs = std::cos(4.0f * float(EIGEN_PI) / size);
  const float sn = std::sin(4.0f * float(EIGEN_PI) / size);
  Matrix3f m;
  m << cs, -sn, 0.0f,
       sn, cs, 0.0f,
       0.0f, 0.0f, 1.0f;
  x[s[1] * 3 + 0] = x[s[0] * 3 + 0];
  x[s[1] * 3 + 1] = -x[s[0] * 3 + 1];
  x[s[1] * 3 + 2] = -x[s[0] * 3 + 2];
  const int hsize = size/2 - 1;
  for (int i = 0; i < hsize; ++i) {
    Eigen::Map<Vector3f>(x.data() + s[(i*2+2) % size] * 3) = m * Eigen::Map<Vector3f>(x.data() + s[i*2] * 3);
  }
  for (int i = hsize; i > 0; --i) {
    Eigen::Map<Vector3f>(x.data() + s[i*2+1] * 3) = m * Eigen::Map<Vector3f>(x.data() + s[(i*2+3) % size] * 3);
  }
}

int Permutation::indexOf(int element) const {
  const auto endIter = std::find(cycle.begin(), cycle.end(), element);
  if (endIter == cycle.end()) { return -1; }
  return int(endIter - cycle.begin());
}

bool Permutation::Normalize(int max_order) {
  order = (int)cycle.size();
  if (order > max_order || max_order % order != 0) {
    std::cout << "ERROR: Invalid cycle length during Normalize." << std::endl;
    return false;
  }
  for (int i = 0; cycle.size() < max_order; ++i) {
    const int e = cycle[i];
    cycle.push_back(e);
  }
  return true;
}

bool Permutation::Apply(VectorXf& x, SymType type) const {
  const std::vector<int>& s = cycle;
  const int fullOrder = (int)cycle.size();
  switch (type) {
  case SymType::Mirror:
    if (order == 1) {
      set_0yz(x, s[0]);
    } else if (order == 2) {
      xref_symmetric(x, s[0], s[1]);
    } else {
      std::cout << "Invalid M symmetry." << std::endl;
      return false;
    }
    break;
  case SymType::Chiral:
    if (order == 1) {
      set_00z(x, s[0]);
    } else if (order == fullOrder) {
      z_symmetric(x, s);
    } else {
      std::cout << "Invalid C symmetry." << std::endl;
      return false;
    }
    break;
  case SymType::Spiegel:
    if (fullOrder % 2 != 0) {
      std::cout << "Invalid S symmetry." << std::endl;
      return false;
    } else if (order == 1) {
      set_000(x, s[0]);
    } else if (order == 2 && fullOrder == 2) {
      pr_symmetric(x, s[0], s[1]);
    } else if (order == 2 && fullOrder > 2) {
      set_00z(x, s[0]);
      pr_symmetric(x, s[0], s[1]);
    } else if (order == fullOrder) {
      z_symmetric(x, s, true);
    } else {
      std::cout << "Invalid S symmetry." << std::endl;
      return false;
    }
    break;
  case SymType::Dihedral:
    if (fullOrder < 4 || fullOrder % 2 != 0) {
      std::cout << "Invalid D symmetry." << std::endl;
      return false;
    } else if (order == 1) {
      set_000(x, s[0]);
    } else if (order == 2 && s[0] == s[2]) {
      set_00z(x, s[0]);
      pr_symmetric(x, s[0], s[1]);
    } else if (order == fullOrder / 2) {
      bool hasMidPoly = false;
      for (int i = 1; i < fullOrder; i += 2) {
        if (s[0] == s[i]) {
          const float angle = (180.0f * (i - 1)) / fullOrder;
          rotZ(x, s[0], angle);
          set_x00(x, s[0]);
          rotZ(x, s[0], -angle);
          z_symmetric(x, s, false, 2);
          hasMidPoly = true;
          break;
        }
      }
      if (!hasMidPoly) {
        std::cout << "Invalid mid-polygon D symmetry." << std::endl;
        return false;
      }
    } else if (order == fullOrder) {
      z_dihedral(x, s);
    } else {
      std::cout << "Invalid D symmetry." << std::endl;
      return false;
    }
    break;
  case SymType::Tetrahedral:
    if (fullOrder != 12) {
      std::cout << "Invalid T symmetry." << std::endl;
      return false;
    } else if (order == 1) {
      set_000(x, s[0]);
    } else if (order == 4) {
      if (s[0] == s[4]) {
        set_xxx(x, s[0]);
        z_dihedral2(x, { s[0], s[1], s[2], s[3] });
      } else if (s[0] == s[5]) {
        set_xnx(x, s[0]);
        z_dihedral2(x, { s[0], s[1], s[2], s[3] });
      } else if (s[0] == s[6]) {
        set_xxn(x, s[0]);
        z_dihedral2(x, { s[0], s[1], s[2], s[3] });
      } else if (s[0] == s[7]) {
        set_xnn(x, s[0]);
        z_dihedral2(x, { s[0], s[1], s[2], s[3] });
      } else {
        std::cout << "Invalid T-4 symmetry." << std::endl;
        return false;
      }
    } else if (order == 6) {
      if (s[0] == s[1]) {
        set_00z(x, s[0]);
        z_offaxis_sym3(x, s[0], s[4], s[8]);
        pr_symmetric(x, s[0], s[2]);
        pr_symmetric(x, s[4], s[5]);
        pr_symmetric(x, s[8], s[9]);
      } else if (s[0] == s[3]) {
        set_0y0(x, s[0]);
        z_offaxis_sym3(x, s[0], s[4], s[8]);
        pr_symmetric(x, s[0], s[1]);
        pr_symmetric(x, s[4], s[5]);
        pr_symmetric(x, s[8], s[10]);
      }  else if (s[0] == s[2]) {
        set_x00(x, s[0]);
        z_offaxis_sym3(x, s[0], s[4], s[8]);
        pr_symmetric(x, s[0], s[1]);
        pr_symmetric(x, s[4], s[6]);
        pr_symmetric(x, s[8], s[9]);
      } else {
        std::cout << "Invalid T-6 symmetry." << std::endl;
        return false;
      }
    } else if (order == 12) {
      z_offaxis_sym3(x, s[0], s[4], s[8]);
      z_dihedral2(x, { s[0], s[1], s[2], s[3] });
      z_dihedral2(x, { s[4], s[5], s[6], s[7] });
      z_dihedral2(x, { s[8], s[9], s[10], s[11] });
    } else {
      std::cout << "Invalid T symmetry." << std::endl;
      return false;
    }
    break;
  case SymType::Octahedral:
    if (fullOrder != 24) {
      std::cout << "Invalid O symmetry." << std::endl;
      return false;
    } else if (order == 1) {
      set_000(x, s[0]);
    } else if (order == 6) {
      if (s[0] == s[1]) {
        set_x00(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        pr_symmetric(x, s[0], s[4]);
        pr_symmetric(x, s[16], s[17]);
        pr_symmetric(x, s[8], s[2]);
      } else if (s[0] == s[5]) {
        set_0y0(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        pr_symmetric(x, s[0], s[1]);
        pr_symmetric(x, s[16], s[20]);
        pr_symmetric(x, s[8], s[9]);
      } else if (s[0] == s[2]) {
        set_00z(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        pr_symmetric(x, s[0], s[1]);
        pr_symmetric(x, s[16], s[17]);
        pr_symmetric(x, s[8], s[12]);
      } else {
        std::cout << "Invalid O-6 symmetry." << std::endl;
        return false;
      }
    } else if (order == 8) {
      if (s[0] == s[8]) {
        set_xxx(x, s[0]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        x_symmetric(x, s[0], s[1]);
        x_symmetric(x, s[2], s[3]);
        x_symmetric(x, s[4], s[5]);
        x_symmetric(x, s[6], s[7]);
      } else if (s[0] == s[12]) {
        set_xnn(x, s[0]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        x_symmetric(x, s[0], s[1]);
        x_symmetric(x, s[2], s[3]);
        x_symmetric(x, s[4], s[5]);
        x_symmetric(x, s[6], s[7]);
      } else if (s[0] == s[13]) {
        set_xxn(x, s[0]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        x_symmetric(x, s[0], s[1]);
        x_symmetric(x, s[2], s[3]);
        x_symmetric(x, s[4], s[5]);
        x_symmetric(x, s[6], s[7]);
      } else if (s[0] == s[9]) {
        set_xnx(x, s[0]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        x_symmetric(x, s[0], s[1]);
        x_symmetric(x, s[2], s[3]);
        x_symmetric(x, s[4], s[5]);
        x_symmetric(x, s[6], s[7]);
      } else {
        std::cout << "Invalid O-8 Symmetry." << std::endl;
        return false;
      }
    } else if (order == 12) {
      if (s[0] == s[3]) {
        set_xx0(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        z_symmetric(x, { s[16], s[22], s[20], s[18] });
        x_symmetric(x, s[16], s[17]);
        x_symmetric(x, s[22], s[23]);
        x_symmetric(x, s[20], s[21]);
        x_symmetric(x, s[18], s[19]);
      } else if (s[0] == s[7]) {
        set_xn0(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        z_symmetric(x, { s[8], s[14], s[12], s[10] });
        z_symmetric(x, { s[16], s[22], s[20], s[18] });
      } else if (s[0] == s[18]) {
        set_x0x(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        z_symmetric(x, { s[8], s[14], s[12], s[10] });
        x_symmetric(x, s[0], s[1]);
        x_symmetric(x, s[2], s[3]);
        x_symmetric(x, s[4], s[5]);
        x_symmetric(x, s[6], s[7]);
      } else if (s[0] == s[23]) {
        set_x0n(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        z_symmetric(x, { s[8], s[14], s[12], s[10] });
        z_symmetric(x, { s[16], s[22], s[20], s[18] });
      } else if (s[0] == s[14]) {
        set_0yy(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        z_symmetric(x, { s[16], s[22], s[20], s[18] });
        x_symmetric(x, s[0], s[1]);
        x_symmetric(x, s[2], s[3]);
        x_symmetric(x, s[4], s[5]);
        x_symmetric(x, s[6], s[7]);
      } else if (s[0] == s[15]) {
        set_0yn(x, s[0]);
        z_offaxis_sym3(x, s[0], s[16], s[8]);
        z_symmetric(x, { s[0], s[6], s[4], s[2] });
        z_symmetric(x, { s[8], s[14], s[12], s[10] });
        z_symmetric(x, { s[16], s[22], s[20], s[18] });
      } else {
        std::cout << "Invalid O-12 symmetry" << std::endl;
        return false;
      }
    } else if (order == 24) {
      z_offaxis_sym3(x, s[0], s[16], s[8]);
      z_symmetric(x, { s[0], s[6], s[4], s[2] });
      z_symmetric(x, { s[8], s[14], s[12], s[10] });
      z_symmetric(x, { s[16], s[22], s[20], s[18] });
      x_symmetric(x, s[0], s[1]);
      x_symmetric(x, s[2], s[3]);
      x_symmetric(x, s[4], s[5]);
      x_symmetric(x, s[6], s[7]);
      x_symmetric(x, s[8], s[9]);
      x_symmetric(x, s[10], s[11]);
      x_symmetric(x, s[12], s[13]);
      x_symmetric(x, s[14], s[15]);
      x_symmetric(x, s[16], s[17]);
      x_symmetric(x, s[18], s[19]);
      x_symmetric(x, s[20], s[21]);
      x_symmetric(x, s[22], s[23]);
    } else {
      std::cout << "Invalid O symmetry." << std::endl;
      return false;
    }
    break;
  case SymType::Icosahedral:
    if (fullOrder != 60) {
      std::cout << "Invalid I symmetry." << std::endl;
      return false;
    } else {
      std::cout << "TODO: Implement I symmetry..." << std::endl;
      return false;
    }
    break;
  case SymType::None:
    return false;
  }
  return true;
}
