#include "symmetry.h"
#include "util.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_set>

void Symmetry::Reset() {
  type = SymType::None;
  order = 0;
  id.clear();
  orders.clear();
  perms.clear();
}

bool Symmetry::Normalize() {
  order = 0;
  for (const Permutation& perm : perms) {
    order = std::max(order, (int)perm.cycle.size());
    orders.insert((int)perm.cycle.size());
  }
  for (Permutation& perm : perms) {
    if (!perm.Normalize(order)) { return false; }
  }
  return true;
}

int Symmetry::Permute(int element, int shift) const {
  for (const Permutation& perm : perms) {
    const int i = perm.indexOf(element);
    if (i >= 0) {
      return perm.cycle[(i + shift) % perm.cycle.size()];
    }
  }
  std::cout << "ERROR: Element doesn't exist in permutation: " << element << std::endl;
  return -1;
}
std::vector<int> Symmetry::Permute(const std::vector<int>& nums, int shift) const {
  std::vector<int> result(nums.begin(), nums.end());
  for (const Permutation& perm : perms) {
    const std::vector<int>& p = perm.cycle;
    for (size_t i = 0; i < p.size(); ++i) {
      const size_t fromIx = p[i];
      const size_t toIx = p[(i + shift) % p.size()];
      result[fromIx] = nums[toIx];
    }
  }
  return result;
}

std::vector<int> Symmetry::IdentityPerm(size_t size) {
  std::vector<int> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = (int)i;
  }
  return result;
}
bool Symmetry::IsIdentityPerm(const std::vector<int>& nums) {
  for (size_t i = 0; i < nums.size(); ++i) {
    if (nums[i] != (int)i) { return false; }
  }
  return true;
}

Symmetry multiply_symmetries(const Symmetry& sym1, const Symmetry& sym2) {
  Symmetry sym_new;
  if (sym2.order != 2) {
    std::cout << "ERROR: Symmetry multiplication must end with an order-2 symmetry!" << std::endl;
    return sym_new;
  }
  std::unordered_set<int> used_ix;
  for (size_t pIx = 0; pIx < sym1.perms.size(); ++pIx) {
    const Permutation& p1 = sym1.perms[pIx];
    if (used_ix.count(p1.rep()) > 0) { continue; }
    Permutation new_perm;
    for (size_t aIx = 0; aIx < p1.cycle.size(); ++aIx) {
      const int a = p1.cycle[aIx];
      new_perm.cycle.push_back(a);
      if (used_ix.count(a) > 0) { return Symmetry(); }
      const int b = sym2.Permute(a);
      new_perm.cycle.push_back(b);
      if (used_ix.count(b) > 0) { return Symmetry(); }
    }
    const size_t prev_size = used_ix.size();
    used_ix.insert(new_perm.cycle.begin(), new_perm.cycle.end());
    new_perm.order = (int)(used_ix.size() - prev_size);
    sym_new.perms.push_back(new_perm);
  }
  sym_new.order = sym1.order * sym2.order;
  return sym_new;
}
Symmetry multiply_symmetries(const Symmetry& sym1, const Symmetry& sym2, const Symmetry& sym3) {
  Symmetry sym_new;
  if (sym3.order != 2) {
    std::cout << "ERROR: Symmetry multiplication must end with an order-2 symmetry!" << std::endl;
    return sym_new;
  }
  std::unordered_set<int> used_ix;
  for (const Permutation& p1 : sym1.perms) {
    if (used_ix.count(p1.rep()) > 0) { continue; }
    Permutation new_perm;
    for (size_t aIx = 0; aIx < sym1.order; ++aIx) {
      int a = p1.cycle[aIx];
      for (size_t bIx = 0; bIx < sym2.order; ++bIx) {
        if (used_ix.count(a) > 0) { return Symmetry(); }
        new_perm.cycle.push_back(a);
        const int b = sym3.Permute(a);
        if (used_ix.count(b) > 0) { return Symmetry(); }
        new_perm.cycle.push_back(b);
        a = sym2.Permute(a);
      }
    }
    const size_t prev_size = used_ix.size();
    used_ix.insert(new_perm.cycle.begin(), new_perm.cycle.end());
    new_perm.order = (int)(used_ix.size() - prev_size);
    sym_new.perms.push_back(new_perm);
  }
  sym_new.order = sym1.order * sym2.order * sym3.order;
  return sym_new;
}
Symmetry multiply_symmetries(const Symmetry& sym1, const Symmetry& sym2, const Symmetry& sym3, const Symmetry& sym4) {
  Symmetry sym_new;
  if (sym3.order != 2 || sym4.order != 2) {
    std::cout << "ERROR: Symmetry multiplication must end with two order-2 symmetries!" << std::endl;
    return sym_new;
  }
  std::unordered_set<int> used_ix;
  for (const Permutation& p1 : sym1.perms) {
    if (used_ix.count(p1.rep()) > 0) { continue; }
    Permutation new_perm;
    for (size_t aIx = 0; aIx < sym1.order; ++aIx) {
      int a = p1.cycle[aIx];
      for (size_t bIx = 0; bIx < sym2.order; ++bIx) {
        if (used_ix.count(a) > 0) { return Symmetry(); }
        new_perm.cycle.push_back(a);
        const int b = sym4.Permute(a);
        if (used_ix.count(b) > 0) { return Symmetry(); }
        new_perm.cycle.push_back(b);
        const int c = sym3.Permute(a);
        if (used_ix.count(c) > 0) { return Symmetry(); }
        new_perm.cycle.push_back(c);
        const int d = sym4.Permute(c);
        if (used_ix.count(d) > 0) { return Symmetry(); }
        new_perm.cycle.push_back(d);
        a = sym2.Permute(a);
      }
    }
    const size_t prev_size = used_ix.size();
    used_ix.insert(new_perm.cycle.begin(), new_perm.cycle.end());
    new_perm.order = (int)(used_ix.size() - prev_size);
    sym_new.perms.push_back(new_perm);
  }
  sym_new.order = sym1.order * sym2.order * sym3.order * sym4.order;
  return sym_new;
}
Symmetry Symmetry::Multiply(const std::vector<Symmetry>& syms) {
  if (syms.size() == 1) {
    return syms[0];
  } else if (syms.size() == 2) {
    return ::multiply_symmetries(syms[0], syms[1]);
  } else if (syms.size() == 3) {
    return ::multiply_symmetries(syms[0], syms[1], syms[2]);
  } else if (syms.size() == 4) {
    return ::multiply_symmetries(syms[0], syms[1], syms[2], syms[3]);
  }
  std::cout << "ERROR: too many symmetries to multiply" << std::endl;
  return Symmetry();
}

Symmetry Symmetry::Dual(const Faces& tris) const {
  if (order == 0) {
    std::cout << "ERROR: Symmetry must be normalized before calling Dual()." << std::endl;
    return Symmetry();
  }
  std::map<int, std::vector<int>> sym_map;
  for (const Permutation& perm : perms) {
    for (int i : perm.cycle) {
      sym_map[i] = perm.cycle;
    }
  }
  std::map<int, Permutation> dual_sym_map;
  for (size_t i = 0; i < tris.size(); ++i) {
    const Face& tri = tris[i];
    const std::vector<int>& v1 = sym_map[tri[0]];
    const std::vector<int>& v2 = sym_map[tri[1]];
    const std::vector<int>& v3 = sym_map[tri[2]];
    const size_t ix1 = std::find(v1.begin(), v1.end(), tri[0]) - v1.begin();
    const size_t ix2 = std::find(v2.begin(), v2.end(), tri[1]) - v2.begin();
    const size_t ix3 = std::find(v3.begin(), v3.end(), tri[2]) - v3.begin();
    int min_u = INT_MAX;
    for (size_t j = 0; j < order; ++j) {
      const int s1 = v1[(ix1 + j) % v1.size()];
      const int s2 = v2[(ix2 + j) % v2.size()];
      const int s3 = v3[(ix3 + j) % v3.size()];
      min_u = std::min(min_u, triangle_id(s1, s2, s3));
    }
    const int maxsize = (int)std::max(v1.size(), std::max(v2.size(), v3.size()));
    if (dual_sym_map.find(min_u) == dual_sym_map.end()) {
      dual_sym_map[min_u] = Permutation({ (int)i });
    } else {
      dual_sym_map[min_u].cycle.push_back((int)i);
      if (dual_sym_map[min_u].cycle.size() > maxsize) {
        std::cout << "WARNING: Something may be wrong with the automorphism..." << std::endl;
      }
    }
  }
  Symmetry dual_sym;
  dual_sym.numElements = (int)tris.size();
  for (const auto& kv : dual_sym_map) {
    dual_sym.perms.push_back(kv.second);
  }
  return dual_sym;
}

bool Symmetry::Load(Symmetry& sym, Symmetry& dsym, const Faces& tris, const char* fname, int line_num, bool zero_indexed) {
  return Load(sym, dsym, tris, fname, std::vector<int>{ line_num }, zero_indexed);
}
bool Symmetry::Load(Symmetry& sym, Symmetry& dsym, const Faces& tris, const char* fname, const std::vector<int>& line_nums, bool zero_indexed) {
  std::vector<Symmetry> syms;
  std::vector<Symmetry> sym_tris;
  for (int line_ix : line_nums) {
    std::ifstream fin(fname);
    if (!fin.is_open()) {
      std::cout << "ERROR: Could not open file: " << fname << std::endl;
      return false;
    }
    std::string line;
    const int zi_sub = (zero_indexed ? 0 : 1);
    int group_line = 0;
    std::map<int, std::vector<int>> sym_map;
    while (true) {
      if (!std::getline(fin, line)) {
        std::cout << "ERROR: Permutation index out of range: " << line_ix << std::endl;
        return false;
      }
      if (line.size() > 0 && line[0] == '#') { continue; }
      group_line += 1;
      if (group_line != line_ix) { continue; }
      ltrim(line);
      rtrim(line);
      std::vector<std::string> strs = split(line, '(');
      std::unordered_set<int> unused_ixs;
      for (size_t fIx = 0; fIx < tris.size(); ++fIx) {
        const Face& face = tris[fIx];
        for (size_t vIx = 0; vIx < face.size(); ++vIx) {
          unused_ixs.insert(face[vIx]);
        }
      }
      size_t max_order = 0;
      for (const std::string& str : strs) {
        if (str.size() == 0) { continue; }
        std::vector<std::string> num_strs = split(str.substr(0, str.size() - 1), ',');
        max_order = std::max(max_order, num_strs.size());
        std::vector<int> nums;
        for (const std::string& str : num_strs) {
          nums.push_back(std::stoi(str) - zi_sub);
        }
        for (int num : nums) {
          sym_map[num] = nums;
          if (unused_ixs.find(num) == unused_ixs.end()) {
            std::cout << "[" << line_ix << "] Found index " << num << " multiple times!" << std::endl;
          } else {
            unused_ixs.erase(num);
          }
        }
      }
      for (auto iter = unused_ixs.begin(); iter != unused_ixs.end(); ++iter) {
        sym_map[*iter].push_back(*iter);
      }
      Symmetry sym;
      for (const auto& kv : sym_map) {
        if (kv.first == kv.second[0]) {
          sym.perms.push_back(Permutation(kv.second));
        }
      }
      if (!sym.Normalize()) { return false; }
      syms.push_back(sym);
      Symmetry sym_tri = sym.Dual(tris);
      if (!sym_tri.Normalize()) { return false; }
      sym_tris.push_back(sym_tri);
      break;
    }
    if (sym_map.empty()) {
      return false;
    }
  }
  sym = Multiply(syms);
  if (sym.Empty()) {
    std::cout << "ERROR: Invalid symmetry multiplication." << std::endl;
    return false;
  }
  dsym = Multiply(sym_tris);
  if (dsym.Empty()) {
    std::cout << "ERROR: Invalid symmetry multiplication." << std::endl;
    return false;
  }
  return true;
}

bool Symmetry::Apply(VectorXf& x) const {
  for (const Permutation& perm : perms) {
    if (!perm.Apply(x, type)) {
      std::cout << "ERROR: Invalid permutation." << std::endl;
      return false;
    }
  }
  return true;
}

std::map<Edge, int> Symmetry::MakeDistMap(const Faces& tris) {
  static const int MAX_IX = 999999;
  static Edges edges;
  static std::map<int, int> i_dist;
  edges.clear();
  make_edges(tris, edges);

  int max_edge = INT_MIN;
  for (const Edge& edge : edges) {
    max_edge = std::max(max_edge, std::max(edge.first, edge.second));
  }
  const int num_verts = max_edge + 1;

  std::map<Edge, int> dist_map;
  for (int i = 0; i < num_verts; ++i) {
    i_dist.clear();
    i_dist[i] = 0;
    while (true) {
      bool changed = false;
      for (const Edge& edge : edges) {
        const int e1 = edge.first;
        const int e2 = edge.second;
        const int p1 = (i_dist.count(e1) > 0 ? i_dist[e1] : MAX_IX);
        const int p2 = (i_dist.count(e2) > 0 ? i_dist[e2] : MAX_IX);
        const int dist = std::min(p1, p2) + 1;
        if (dist < p1) {
          i_dist[e1] = dist;
          changed = true;
        }
        if (dist < p2) {
          i_dist[e2] = dist;
          changed = true;
        }
      }
      if (!changed) { break; }
    }
    for (const auto& kv : i_dist) {
      dist_map[Edge(i, kv.first)] = kv.second;
    }
  }

  return dist_map;
}

std::string Symmetry::MakeId(const Faces& tris, const std::map<Edge, int>& dist_map) {
  static std::unordered_set<int> unused_ixs;
  static std::map<std::pair<int, int>, int> swap_dist_count;

  unused_ixs.clear();
  for (size_t fIx = 0; fIx < tris.size(); ++fIx) {
    const Face& face = tris[fIx];
    for (size_t vIx = 0; vIx < face.size(); ++vIx) {
      unused_ixs.insert(face[vIx]);
    }
  }

  order = 0;
  swap_dist_count.clear();
  for (const Permutation& perm : perms) {
    order = std::max(order, (int)perm.cycle.size());
    for (size_t i = 0; i < perm.cycle.size(); ++i) {
      const int fromIx = perm.cycle[i];
      const int toIx = perm.cycle[(i + 1) % perm.cycle.size()];
      const int dist = dist_map.at(Edge(fromIx, toIx));
      swap_dist_count[std::pair((int)perm.cycle.size(), dist)] += 1;
      unused_ixs.erase(fromIx);
    }
  }
  std::stringstream ss;
  ss << std::setw(2) << std::setfill(' ') << order << "(";
  if (unused_ixs.size() > 0) {
    ss << "1a:" << unused_ixs.size() << "  ";
  }
  bool first = true;
  for (const auto& kv : swap_dist_count) {
    ss << kv.first.first << char('a' + kv.first.second) << ":" << kv.second << " ";
  }
  id = ss.str();
  id = id.substr(0, id.size()-1) + ")";
  return id;
}

bool has_sym(const Verts3D& x, const Symmetry& sym, Vector3f axis, bool alternate) {
  axis.normalize();
  for (const Permutation& perm : sym.perms) {
    const std::vector<int>& p = perm.cycle;
    Vector3f sum = Vector3f::Zero();
    Vector3f same_cp = Vector3f::Zero();
    float same_mag = 0.0f;
    bool first = true;
    for (size_t i = 0; i < p.size(); ++i) {
      Vector3f a = x[p[i]];
      Vector3f b = x[p[(i + 1) % p.size()]];
      Vector3f c = x[p[(i + 2) % p.size()]];
      if (alternate) {
        if (i % 2 == 1) {
          a -= axis * (2.0f * a.dot(axis));
          c -= axis * (2.0f * c.dot(axis));
        } else {
          b -= axis * (2.0f * b.dot(axis));
        }
      }
      const Vector3f cp = (b - a).cross(c - b);
      if (std::abs((b - a).dot(axis)) > 1e-3f) {
        return false;
      }
      const float mag = (b - a).norm();
      sum += a;
      if (first) {
        same_cp = cp;
        same_mag = mag;
        first = false;
      } else if ((cp - same_cp).norm() > 1e-3f || std::abs(mag - same_mag) > 1e-3f) {
        return false;
      }
    }
    const float cp = sum.cross(axis).norm();
    if (cp > 1e-3f) {
      return false;
    }
  }
  return true;
}
std::string Symmetry::Test(const Verts3D& x) const {
  static const char axisNames[] = "xyzd";
  static const Vector3f testAxes[] = {
    Vector3f::UnitX(),
    Vector3f::UnitY(),
    Vector3f::UnitZ(),
    Vector3f::Ones(),
  };
  for (int i = 0; i < 4; ++i) {
    if (has_sym(x, *this, testAxes[i], false)) {
      return axisNames[i] + std::string("C") + std::to_string(order);
    } else if (order % 2 == 0 && has_sym(x, *this, testAxes[i], true)) {
      return axisNames[i] + std::string("S") + std::to_string(order);
    }
  }
  return "";
}

SymType Symmetry::GetType(char c) {
  switch (c) {
  case 'M': return SymType::Mirror;
  case 'C': return SymType::Chiral;
  case 'S': return SymType::Spiegel;
  case 'D': return SymType::Dihedral;
  case 'T': return SymType::Tetrahedral;
  case 'O': return SymType::Octahedral;
  case 'I': return SymType::Icosahedral;
  default: return SymType::None;
  }
}

int Symmetry::GetProdSize(SymType symType) {
  switch (symType) {
  case SymType::Mirror: return 1;
  case SymType::Chiral: return 1;
  case SymType::Spiegel: return 1;
  case SymType::Dihedral: return 2;
  case SymType::Tetrahedral: return 3;
  case SymType::Octahedral: return 3;
  case SymType::Icosahedral: return 4;
  default: return 0;
  }
}

void SymmetryList::Reset() {
  syms.clear();
  symDuals.clear();
  symsByOrder.clear();
  symsById.clear();
}

bool SymmetryList::Load(const Faces& tris, const char* fname, bool zero_indexed) {
  Reset();
  std::cout << "Creating Distance map..." << std::endl;
  const std::map<Edge, int> dist_map = Symmetry::MakeDistMap(tris);
  std::cout << "Loading all automorphisms from " << fname << "..." << std::endl;
  std::ifstream fin(fname);
  std::string line;
  const int zi_sub = (zero_indexed ? 0 : 1);
  int group_line = 0;
  while (std::getline(fin, line)) {
    if (line.size() > 0 && line[0] == '#') { continue; }
    group_line += 1;
    ltrim(line);
    rtrim(line);
    std::map<int, std::vector<int>> sym_map;
    std::vector<std::string> strs = split(line, '(');
    std::unordered_set<int> unused_ixs;
    for (size_t fIx = 0; fIx < tris.size(); ++fIx) {
      const Face& face = tris[fIx];
      for (size_t vIx = 0; vIx < face.size(); ++vIx) {
        unused_ixs.insert(face[vIx]);
      }
    }
    const int numElements = (int)unused_ixs.size();
    size_t max_order = 0;
    for (const std::string& str : strs) {
      if (str.size() == 0) { continue; }
      std::vector<std::string> num_strs = split(str.substr(0, str.size() - 1), ',');
      max_order = std::max(max_order, num_strs.size());
      std::vector<int> nums;
      for (const std::string& str : num_strs) {
        nums.push_back(std::stoi(str) - zi_sub);
      }
      for (int num : nums) {
        sym_map[num] = nums;
        if (unused_ixs.find(num) == unused_ixs.end()) {
          std::cout << "ERROR: Found index " << num << " multiple times!" << std::endl;
          return false;
        } else {
          unused_ixs.erase(num);
        }
      }
    }
    for (auto iter = unused_ixs.begin(); iter != unused_ixs.end(); ++iter) {
      sym_map[*iter].push_back(*iter);
    }
    Symmetry sym;
    sym.numElements = numElements;
    for (const auto& kv : sym_map) {
      if (kv.first == kv.second[0]) {
        sym.perms.push_back(Permutation(kv.second));
      }
    }
    symsById[sym.MakeId(tris, dist_map)].push_back((int)syms.size());
    if (!sym.Normalize()) { return false; }
    Symmetry sym_tri = sym.Dual(tris);
    if (!sym_tri.Normalize()) { return false; }
    symsByOrder[sym.order].push_back((int)syms.size());
    syms.push_back(sym);
    symDuals.push_back(sym_tri);
  }
  std::cout << "Found " << syms.size() << " non-trivial automorphisms with:" << std::endl;
  std::cout << "  " << symsByOrder.size() << " different orders and" << std::endl;
  std::cout << "  " << symsById.size() << " unique permutation types" << std::endl;
  return true;
}

bool IsChiralPossible(const Symmetry& rep) {
  //Check if chiral can be supported in the most basic form
  if (rep.orders.size() != 1 && (rep.orders.size() != 2 || rep.orders.count(1) == 0)) {
    return false;
  }
  return true;
}
bool IsChiralEmbeddable(const Symmetry& rep, const Faces& tris, const std::set<Edge>& edgeSet) {
  //Check for crossed edges on a cycle
  for (const Permutation& perm : rep.perms) {
    int numEdges = 0;
    for (int i = 0; i < perm.order; ++i) {
      for (int j = 2; j <= perm.order - 2; ++j) {
        const int p1 = perm.cycle[i];
        const int p2 = perm.cycle[(i + j + perm.order) % perm.order];
        if (edgeSet.count(Edge(p1, p2)) > 0 || edgeSet.count(Edge(p2, p1)) > 0) {
          numEdges += 1;
        }
      }
    }
    if (numEdges > 2) {
      return false;
    }
  }
  //Check for degenerate triangles
  for (const Face& tri : tris) {
    const int p0 = rep.Permute(tri[0]);
    const int p1 = rep.Permute(tri[1]);
    const int p2 = rep.Permute(tri[2]);
    if ((tri[0] == p0 && (tri[1] == p1 || tri[2] == p2 || (tri[1] == p2 && tri[2] == p1))) ||
        (tri[1] == p1 && (tri[0] == p0 || tri[2] == p2 || (tri[2] == p0 && tri[0] == p2))) ||
        (tri[2] == p2 && (tri[0] == p0 || tri[1] == p1 || (tri[0] == p1 && tri[1] == p0)))) {
      return false;
    }
  }
  return true;
}
bool GetChiralEmbeddableIx(const std::vector<int>& ixs, const std::vector<Symmetry>& syms, const Faces& tris, const std::set<Edge>& edgeSet, int& repIx) {
  for (int symIx : ixs) {
    if (IsChiralEmbeddable(syms[symIx], tris, edgeSet)) {
      repIx = symIx;
      return true;
    }
  }
  return false;
}

bool IsMirrorEmbeddable(const Symmetry& rep, const Faces& tris, const std::set<Edge>& edgeSet) {
  //Check for degenerate triangles
  for (const Face& tri : tris) {
    const int p0 = rep.Permute(tri[0]);
    const int p1 = rep.Permute(tri[1]);
    const int p2 = rep.Permute(tri[2]);
    //Check if triangle is entirely on the mirror plane
    if (tri[0] == p0 && tri[1] == p1 && tri[2] == p2) { return false; }
    //Check if triangle spans the gap
    if (tri[0] == p1 && tri[1] == p0 && tri[2] != p2) { return false; }
    if (tri[0] == p2 && tri[2] == p0 && tri[1] != p1) { return false; }
    if (tri[1] == p2 && tri[2] == p1 && tri[0] != p0) { return false; }
  }
  //Symmetry must be splittable into two halves
  int startIx = 0;
  while (startIx == rep.Permute(startIx)) { ++startIx; }
  std::vector<int> to_explore;
  std::unordered_set<int> explored;
  to_explore.push_back(startIx);
  while (!to_explore.empty()) {
    const int ix = to_explore.back();
    to_explore.pop_back();
    if (explored.count(ix) > 0) { continue; }
    explored.insert(ix);
    for (const Face& tri : tris) {
      int a, b;
      if (tri[0] == ix) {
        a = tri[1]; b= tri[2];
      } else if (tri[1] == ix) {
        a = tri[0]; b = tri[2];
      } else if (tri[2] == ix) {
        a = tri[0]; b = tri[1];
      } else {
        continue;
      }
      const int pa = rep.Permute(a);
      const int pb = rep.Permute(b);
      if (explored.count(a) == 0 && pa != a && pa != ix) {
        if (explored.count(pa) > 0) { return false; }
        to_explore.push_back(a);
      }
      if (explored.count(b) == 0 && pb != b && pb != ix) {
        if (explored.count(pb) > 0) { return false; }
        to_explore.push_back(b);
      }
    }
  }
  return true;
}
bool GetMirrorEmbeddableIx(const std::vector<int>& ixs, const std::vector<Symmetry>& syms, const Faces& tris, const std::set<Edge>& edgeSet, int& repIx) {
  for (int symIx : ixs) {
    if (IsMirrorEmbeddable(syms[symIx], tris, edgeSet)) {
      repIx = symIx;
      return true;
    }
  }
  return false;
}

bool IsSpiegelPossible(const Symmetry& rep) {
  //Check if Spiegel can be supported in the most basic form
  if (rep.order % 2 != 0) { return false; }
  if (rep.orders.size() > 2) { return false; }
  if (rep.orders.count(1) > 0) { return false; }
  if (rep.orders.size() == 2 && rep.orders.count(2) == 0) { return false; }
  return true;
}
bool IsSpiegelEmbeddable(const Symmetry& rep, const Faces& tris, const std::set<Edge>& edgeSet) {
  //Check for crossed edges on a cycle
  for (const Permutation& perm : rep.perms) {
    int numEdges = 0;
    for (int i = 0; i < perm.order/2; ++i) {
      for (int j = 2; j <= perm.order/2 - 2; ++j) {
        const int p1 = perm.cycle[i*2];
        const int p2 = perm.cycle[((i + j)*2 + perm.order) % perm.order];
        if (edgeSet.count(Edge(p1, p2)) > 0 || edgeSet.count(Edge(p2, p1)) > 0) {
          numEdges += 1;
        }
      }
    }
    if (numEdges > 2) {
      return false;
    }
  }
  //Check if any edges cross through the center
  for (const Permutation& perm : rep.perms) {
    const int halfOrder = perm.order / 2;
    if (halfOrder % 2 == 1) {
      for (int i = 0; i < halfOrder; ++i) {
        const int p1 = perm.cycle[i];
        const int p2 = perm.cycle[i + halfOrder];
        if (edgeSet.count(Edge(p1, p2)) > 0 || edgeSet.count(Edge(p2, p1)) > 0) {
          return false;
        }
      }
    }
  }
  return true;
}
bool GetSpiegelEmbeddableIx(const std::vector<int>& ixs, const std::vector<Symmetry>& syms, const Faces& tris, const std::set<Edge>& edgeSet, int& repIx) {
  for (int symIx : ixs) {
    if (IsSpiegelEmbeddable(syms[symIx], tris, edgeSet)) {
      repIx = symIx;
      return true;
    }
  }
  return false;
}

void SymmetryList::Analyze(const Faces& tris, bool embeddable, std::ostream& os) const {
  //Preprocessing
  Edges edges;
  make_edges(tris, edges);
  std::set<Edge> edgeSet(edges.begin(), edges.end());
  //FaceMap faceMap;
  //for (const Face& tri : tris) {
  //  faceMap[triangle_id(tri)] = tri;
  //}

  //Find unique M symmetries
  for (const auto& kv : symsById) {
    int repIx = kv.second[0];
    const Symmetry& rep = syms[repIx];
    if (rep.order != 2) { continue; }
    if (!IsChiralPossible(rep)) { continue; }
    const bool repE = GetMirrorEmbeddableIx(kv.second, syms, tris, edgeSet, repIx);
    if (!repE && embeddable) { continue; }
    os << (repE ? " " : "*");
    os << "M" << rep.order << " : " << (repIx + 1) << "   " << rep.id << std::endl;
  }

  //Find unique C symmetries
  for (const auto& kv : symsById) {
    int repIx = kv.second[0];
    const Symmetry& rep = syms[repIx];
    if (!IsChiralPossible(rep)) { continue; }
    const bool repE = GetChiralEmbeddableIx(kv.second, syms, tris, edgeSet, repIx);
    if (!repE && embeddable) { continue; }
    os << (repE ? " " : "*");
    os << "C" << rep.order << " : " << (repIx + 1) << "   " << rep.id << std::endl;
  }

  //Find unique S symmetries
  for (const auto& kv : symsById) {
    int repIx = kv.second[0];
    const Symmetry& rep = syms[repIx];
    if (!IsSpiegelPossible(rep)) { continue; }
    const bool repE = GetSpiegelEmbeddableIx(kv.second, syms, tris, edgeSet, repIx);
    if (!repE && embeddable) { continue; }
    os << (repE ? " " : "*");
    os << "S" << rep.order << " : " << (repIx + 1) << "   " << rep.id << std::endl;
  }

  //Find unique D symmetries
  for (const auto& kv : symsById) {
    int repIx = kv.second[0];
    if (!IsChiralPossible(syms[repIx])) { continue; }
    const bool repE = GetChiralEmbeddableIx(kv.second, syms, tris, edgeSet, repIx);
    if (!repE && embeddable) { continue; }
    const Symmetry& rep = syms[repIx];
    const Symmetry& rep_d = symDuals[repIx];

    for (const auto& kv2 : symsById) {
      int repIx2 = kv2.second[0];
      const Symmetry& rep2 = syms[repIx2];
      if (rep2.order != 2) { continue; }
      if (rep.order == 2 && rep.id < rep2.id) { continue; }
      if (!IsChiralPossible(rep2)) { continue; }
      const bool rep2E = GetChiralEmbeddableIx(kv.second, syms, tris, edgeSet, repIx2);
      if (!rep2E && embeddable) { continue; }

      for (const int symIx : kv2.second) {
        const Symmetry& sym = syms[symIx];
        const Symmetry& sym_d = symDuals[symIx];
        if (embeddable && !IsChiralEmbeddable(sym, tris, edgeSet)) { continue; }

        //Look for D cycles
        std::vector<int> perm = Symmetry::IdentityPerm(sym.numElements);
        perm = rep.Permute(perm);
        perm = sym.Permute(perm);
        if (Symmetry::IsIdentityPerm(perm)) { continue; }
        perm = rep.Permute(perm);
        perm = sym.Permute(perm);
        if (!Symmetry::IsIdentityPerm(perm)) { continue; }

        //Try multiplying the symmetries to see if the product is valid
        Symmetry prod = Symmetry::Multiply({rep, sym});
        if (prod.Empty()) { continue; }
        Symmetry dualProd = Symmetry::Multiply({ rep_d, sym_d });
        if (dualProd.Empty()) { continue; }

        //Check if the product has a valid order
        bool hasFullOrder = false;
        for (const Permutation& prodPerm : prod.perms) {
          if (prodPerm.order == prod.order) {
            hasFullOrder = true;
            break;
          }
        }
        if (!hasFullOrder) { continue; }

        //Also check to see if a poles make sense
        bool hasValidPole = true;
        if (prod.order >= 6) {
          for (const Permutation& prodPerm : prod.perms) {
            if (prodPerm.order == prod.order/2) {
              if (prodPerm.cycle[0] == prodPerm.cycle[2] || prodPerm.cycle[0] == prodPerm.cycle[4]) {
                hasValidPole = false;
                break;
              }
            } else if (prodPerm.order == 2) {
              if (prodPerm.cycle[0] != prodPerm.cycle[2] || prodPerm.cycle[0] != prodPerm.cycle[4]) {
                hasValidPole = false;
                break;
              }
            } else if (prodPerm.order != prod.order){
              hasValidPole = false;
              break;
            }
          }
        }
        if (!hasValidPole) { continue; }

        //Print result
        os << (repE && rep2E ? " " : "*");
        os << "D" << (rep.order) << " : " << (repIx + 1) << "," << (symIx + 1)
                  << "   " << rep.id << " x" << sym.id << std::endl;
        break;
      }
    }
  }

  //Find unique T symmetries
  for (const auto& kv3 : symsById) {
    int rep3Ix = kv3.second[0];
    if (syms[rep3Ix].order != 3) { continue; }
    if (!IsChiralPossible(syms[rep3Ix])) { continue; }
    const bool rep3E = GetChiralEmbeddableIx(kv3.second, syms, tris, edgeSet, rep3Ix);
    if (!rep3E && embeddable) { continue; }
    const Symmetry& rep3 = syms[rep3Ix];
    const Symmetry& rep3_d = symDuals[rep3Ix];

    for (const auto& kv2a : symsById) {
      int repIx2a = kv2a.second[0];
      const Symmetry& rep2a = syms[repIx2a];
      if (rep2a.order != 2) { continue; }
      if (!IsChiralPossible(rep2a)) { continue; }
      const bool rep2aE = GetChiralEmbeddableIx(kv2a.second, syms, tris, edgeSet, repIx2a);
      if (!rep2aE && embeddable) { continue; }

      for (const auto& kv2b : symsById) {
        int repIx2b = kv2b.second[0];
        const Symmetry& rep2b = syms[repIx2b];
        if (rep2b.order != 2) { continue; }
        if (rep2a.id < rep2b.id) { continue; }
        if (!IsChiralPossible(rep2b)) { continue; }
        const bool rep2bE = GetChiralEmbeddableIx(kv2b.second, syms, tris, edgeSet, repIx2b);
        if (!rep2bE && embeddable) { continue; }

        //At this point we have a compatible set of symmetries, find any that work out
        bool found = false;
        for (const int symAIx : kv2a.second) {
          const Symmetry& symA = syms[symAIx];
          const Symmetry& symA_d = symDuals[symAIx];
          if (embeddable && !IsChiralEmbeddable(symA, tris, edgeSet)) { continue; }
          std::vector<int> perm = Symmetry::IdentityPerm(symA.numElements);

          //3-A
          perm = symA.Permute(perm);
          perm = rep3.Permute(perm);
          if (Symmetry::IsIdentityPerm(perm)) { continue; }
          perm = symA.Permute(perm);
          perm = rep3.Permute(perm);
          if (Symmetry::IsIdentityPerm(perm)) { continue; }
          perm = symA.Permute(perm);
          perm = rep3.Permute(perm);
          if (!Symmetry::IsIdentityPerm(perm)) { continue; }

          for (const int symBIx : kv2b.second) {
            const Symmetry& symB = syms[symBIx];
            const Symmetry& symB_d = symDuals[symBIx];
            if (embeddable && !IsChiralEmbeddable(symB, tris, edgeSet)) { continue; }
            perm = Symmetry::IdentityPerm(symB.numElements);
            //3-B
            perm = symB.Permute(perm);
            perm = rep3.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = symB.Permute(perm);
            perm = rep3.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = symB.Permute(perm);
            perm = rep3.Permute(perm);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }

            //A-B
            perm = symA.Permute(perm);
            perm = symB.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = symA.Permute(perm);
            perm = symB.Permute(perm);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }

            //Directionality
            perm = rep3.Permute(perm);
            perm = symA.Permute(perm);
            perm = symB.Permute(perm);
            perm = rep3.Permute(perm,2);
            perm = symA.Permute(perm);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }

            //Try multiplying the symmetries to see if the product is valid
            Symmetry prod = Symmetry::Multiply({ rep3, symA, symB });
            if (prod.Empty()) { continue; }
            Symmetry dualProd = Symmetry::Multiply({ rep3_d, symA_d, symB_d });
            if (dualProd.Empty()) { continue; }

            //Check if the product has a valid order
            bool hasValidOrder = true;
            bool hasFullOrder = (prod.perms.size() == 1);
            for (const Permutation& prodPerm : prod.perms) {
              if (prodPerm.order == 12) { hasFullOrder = true; }
              if (prodPerm.order != 12 && prodPerm.order != 6 && prodPerm.order != 4 && prodPerm.order != 1) {
                hasValidOrder = false;
                break;
              }
            }
            if (!hasValidOrder) { continue; }
            if (embeddable && !hasFullOrder) { continue; }

            //Print result
            os << (rep3E && rep2aE && rep2bE && hasFullOrder ? " " : "*");
            os << "T  : " << (rep3Ix + 1) << "," << (symAIx + 1) << "," << (symBIx + 1)
              << "   " << rep3.id << " x" << symA.id << " x" << symB.id << std::endl;
            found = true;
            break;
          }
          if (found) { break; }
        }
      }
    }
  }

  //Find unique O symmetries
  for (const auto& kv4 : symsById) {
    int rep4Ix = kv4.second[0];
    if (syms[rep4Ix].order != 4) { continue; }
    if (!IsChiralPossible(syms[rep4Ix])) { continue; }
    const bool rep4E = GetChiralEmbeddableIx(kv4.second, syms, tris, edgeSet, rep4Ix);
    if (!rep4E && embeddable) { continue; }
    const Symmetry& rep4 = syms[rep4Ix];
    const Symmetry& rep4_d = symDuals[rep4Ix];

    for (const auto& kv3 : symsById) {
      int repIx3 = kv3.second[0];
      if (syms[repIx3].order != 3) { continue; }
      if (!IsChiralPossible(syms[repIx3])) { continue; }
      const bool rep3E = GetChiralEmbeddableIx(kv3.second, syms, tris, edgeSet, repIx3);
      if (!rep3E && embeddable) { continue; }

      for (const auto& kv2 : symsById) {
        int repIx2 = kv2.second[0];
        if (syms[repIx2].order != 2) { continue; }
        if (!IsChiralPossible(syms[repIx2])) { continue; }
        const bool rep2E = GetChiralEmbeddableIx(kv2.second, syms, tris, edgeSet, repIx2);
        if (!rep2E && embeddable) { continue; }

        //At this point we have a compatible set of symmetries, find any that work out
        bool found = false;
        for (const int sym3Ix : kv3.second) {
          const Symmetry& sym3 = syms[sym3Ix];
          const Symmetry& sym3_d = symDuals[sym3Ix];
          if (embeddable && !IsChiralEmbeddable(sym3, tris, edgeSet)) { continue; }
          std::vector<int> perm = Symmetry::IdentityPerm(sym3.numElements);

          //O
          perm = rep4.Permute(perm);
          perm = sym3.Permute(perm);
          if (Symmetry::IsIdentityPerm(perm)) { continue; }
          perm = rep4.Permute(perm);
          perm = sym3.Permute(perm);
          if (Symmetry::IsIdentityPerm(perm)) { continue; }
          perm = rep4.Permute(perm);
          perm = sym3.Permute(perm);
          if (Symmetry::IsIdentityPerm(perm)) { continue; }
          perm = rep4.Permute(perm);
          perm = sym3.Permute(perm);
          if (!Symmetry::IsIdentityPerm(perm)) { continue; }

          for (const int sym2Ix : kv2.second) {
            const Symmetry& sym2 = syms[sym2Ix];
            const Symmetry& sym2_d = symDuals[sym2Ix];
            if (embeddable && !IsChiralEmbeddable(sym2, tris, edgeSet)) { continue; }
            perm = Symmetry::IdentityPerm(sym2.numElements);

            //D4
            perm = rep4.Permute(perm);
            perm = sym2.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = rep4.Permute(perm);
            perm = sym2.Permute(perm);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }
            
            //T
            perm = sym3.Permute(perm);
            perm = sym2.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = sym3.Permute(perm);
            perm = sym2.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = sym3.Permute(perm);
            perm = sym2.Permute(perm);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }

            //Directionality
            perm = sym3.Permute(perm);
            perm = rep4.Permute(perm);
            perm = sym2.Permute(perm);
            perm = sym3.Permute(perm);
            perm = sym2.Permute(perm);
            perm = rep4.Permute(perm,3);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }

            //Try multiplying the symmetries to see if the product is valid
            Symmetry prod = Symmetry::Multiply({ sym3, rep4, sym2 });
            if (prod.Empty()) { continue; }
            Symmetry dualProd = Symmetry::Multiply({ sym3_d, rep4_d, sym2_d });
            if (dualProd.Empty()) { continue; }

            //Check if the product has a valid order
            bool hasValidOrder = true;
            bool hasFullOrder = (prod.perms.size() == 1);
            for (const Permutation& prodPerm : prod.perms) {
              if (prodPerm.order == 24) { hasFullOrder = true; }
              if (prodPerm.order != 24 && prodPerm.order != 12 && prodPerm.order != 8 && prodPerm.order != 6 && prodPerm.order != 1) {
                hasValidOrder = false;
                break;
              }
            }
            if (!hasValidOrder) { continue; }
            if (embeddable && !hasFullOrder) { continue; }

            //Also make sure there are no triangles going through the center
            bool trianglesE = true;
            if (embeddable) {
              for (const Permutation& prodPerm : prod.perms) {
                if (prodPerm.order == 12) {
                  const int maxV = *std::max_element(prodPerm.cycle.begin(), prodPerm.cycle.end());
                  VectorXf x((maxV+1) * 3);
                  for (int i : prodPerm.cycle) {
                    x[i*3+0] = 1;
                    x[i*3+1] = 2;
                    x[i*3+2] = 3;
                  }
                  prodPerm.Apply(x, SymType::Octahedral);
                  const std::unordered_set permSet(prodPerm.cycle.begin(), prodPerm.cycle.end());
                  for (const Face& tri : tris) {
                    if (permSet.count(tri[0]) > 0 && permSet.count(tri[1]) > 0 && permSet.count(tri[2]) > 0) {
                      const float centerX = x[tri[0]*3+0] + x[tri[1]*3+0] + x[tri[2]*3+0];
                      const float centerY = x[tri[0]*3+1] + x[tri[1]*3+1] + x[tri[2]*3+1];
                      const float centerZ = x[tri[0]*3+2] + x[tri[1]*3+2] + x[tri[2]*3+2];
                      if (centerX * centerX + centerY * centerY + centerZ * centerZ < 1e-3f) {
                        trianglesE = false;
                      }
                    }
                  }
                  if (!trianglesE) { break; }
                }
              }
              if (!trianglesE) { continue; }
            }

            //Print result
            os << (rep4E && rep3E && rep2E && trianglesE && hasFullOrder ? " " : "*");
            os << "O  : " << (sym3Ix + 1) << "," << (rep4Ix + 1) << "," << (sym2Ix + 1)
                      << "   " << sym3.id << " x" << rep4.id << " x" << sym2.id << std::endl;
            found = true;
            break;
          }
          if (found) { break; }
        }
      }
    }
  }

  //Find unique I symmetries
  for (const auto& kv5 : symsById) {
    int rep5Ix = kv5.second[0];
    if (syms[rep5Ix].order != 5) { continue; }
    if (!IsChiralPossible(syms[rep5Ix])) { continue; }
    const bool rep5E = GetChiralEmbeddableIx(kv5.second, syms, tris, edgeSet, rep5Ix);
    if (!rep5E && embeddable) { continue; }
    const Symmetry& rep5 = syms[rep5Ix];
    const Symmetry& rep5_d = symDuals[rep5Ix];

    for (const auto& kv3 : symsById) {
      int repIx3 = kv3.second[0];
      if (syms[repIx3].order != 3) { continue; }
      if (!IsChiralPossible(syms[repIx3])) { continue; }
      const bool rep3E = GetChiralEmbeddableIx(kv3.second, syms, tris, edgeSet, repIx3);
      if (!rep3E && embeddable) { continue; }

      for (const auto& kv2z : symsById) {
        int repIx2z = kv2z.second[0];
        if (syms[repIx2z].order != 2) { continue; }
        if (!IsChiralPossible(syms[repIx2z])) { continue; }
        const bool rep2zE = GetChiralEmbeddableIx(kv2z.second, syms, tris, edgeSet, repIx2z);
        if (!rep2zE && embeddable) { continue; }

        for (const auto& kv2x : symsById) {
          int repIx2x = kv2x.second[0];
          if (syms[repIx2x].order != 2) { continue; }
          if (!IsChiralPossible(syms[repIx2x])) { continue; }
          const bool rep2xE = GetChiralEmbeddableIx(kv2x.second, syms, tris, edgeSet, repIx2x);
          if (!rep2xE && embeddable) { continue; }

          //At this point we have a compatible set of symmetries, find any that work out
          bool found = false;
          for (const int sym3Ix : kv3.second) {
            const Symmetry& sym3 = syms[sym3Ix];
            const Symmetry& sym3_d = symDuals[sym3Ix];
            if (embeddable && !IsChiralEmbeddable(sym3, tris, edgeSet)) { continue; }
            std::vector<int> perm = Symmetry::IdentityPerm(sym3.numElements);

            //I
            perm = rep5.Permute(perm);
            perm = sym3.Permute(perm);
            if (Symmetry::IsIdentityPerm(perm)) { continue; }
            perm = rep5.Permute(perm);
            perm = sym3.Permute(perm);
            if (!Symmetry::IsIdentityPerm(perm)) { continue; }

            for (const int sym2zIx : kv2z.second) {
              const Symmetry& sym2z = syms[sym2zIx];
              const Symmetry& sym2z_d = symDuals[sym2zIx];
              if (embeddable && !IsChiralEmbeddable(sym2z, tris, edgeSet)) { continue; }
              perm = Symmetry::IdentityPerm(sym2z.numElements);

              //D5
              perm = rep5.Permute(perm);
              perm = sym2z.Permute(perm);
              if (Symmetry::IsIdentityPerm(perm)) { continue; }
              perm = rep5.Permute(perm);
              perm = sym2z.Permute(perm);
              if (Symmetry::IsIdentityPerm(perm)) { continue; }
              perm = rep5.Permute(perm);
              perm = sym2z.Permute(perm);
              if (Symmetry::IsIdentityPerm(perm)) { continue; }
              perm = rep5.Permute(perm);
              perm = sym2z.Permute(perm);
              if (Symmetry::IsIdentityPerm(perm)) { continue; }
              perm = rep5.Permute(perm);
              perm = sym2z.Permute(perm);
              if (!Symmetry::IsIdentityPerm(perm)) { continue; }

              //D3
              perm = sym3.Permute(perm);
              perm = sym2z.Permute(perm);
              if (Symmetry::IsIdentityPerm(perm)) { continue; }
              perm = sym3.Permute(perm);
              perm = sym2z.Permute(perm);
              if (Symmetry::IsIdentityPerm(perm)) { continue; }
              perm = sym3.Permute(perm);
              perm = sym2z.Permute(perm);
              if (!Symmetry::IsIdentityPerm(perm)) { continue; }

              for (const int sym2xIx : kv2x.second) {
                const Symmetry& sym2x = syms[sym2xIx];
                const Symmetry& sym2x_d = symDuals[sym2xIx];
                if (embeddable && !IsChiralEmbeddable(sym2x, tris, edgeSet)) { continue; }
                perm = Symmetry::IdentityPerm(sym2x.numElements);

                //T
                perm = rep5.Permute(perm);
                perm = sym2x.Permute(perm);
                if (Symmetry::IsIdentityPerm(perm)) { continue; }
                perm = rep5.Permute(perm);
                perm = sym2x.Permute(perm);
                if (Symmetry::IsIdentityPerm(perm)) { continue; }
                perm = rep5.Permute(perm);
                perm = sym2x.Permute(perm);
                if (!Symmetry::IsIdentityPerm(perm)) { continue; }

                //T
                perm = sym3.Permute(perm);
                perm = sym2x.Permute(perm);
                if (Symmetry::IsIdentityPerm(perm)) { continue; }
                perm = sym3.Permute(perm);
                perm = sym2x.Permute(perm);
                if (Symmetry::IsIdentityPerm(perm)) { continue; }
                perm = sym3.Permute(perm);
                perm = sym2x.Permute(perm);
                if (!Symmetry::IsIdentityPerm(perm)) { continue; }

                //D2
                perm = sym2z.Permute(perm);
                perm = sym2x.Permute(perm);
                if (Symmetry::IsIdentityPerm(perm)) { continue; }
                perm = sym2z.Permute(perm);
                perm = sym2x.Permute(perm);
                if (!Symmetry::IsIdentityPerm(perm)) { continue; }

                //Directionality
                perm = rep5.Permute(perm);
                perm = sym3.Permute(perm);
                perm = sym2z.Permute(perm);
                perm = sym2x.Permute(perm);
                perm = rep5.Permute(perm, 4);
                perm = sym3.Permute(perm);
                perm = sym2x.Permute(perm);
                if (!Symmetry::IsIdentityPerm(perm)) { continue; }

                //Try multiplying the symmetries to see if the product actually works
                Symmetry prod = Symmetry::Multiply({ sym3, rep5, sym2z, sym2x });
                if (prod.Empty()) { continue; }
                Symmetry dualProd = Symmetry::Multiply({ sym3_d, rep5_d, sym2z_d, sym2x_d });
                if (dualProd.Empty()) { continue; }

                //Check if the product has a valid order
                bool hasValidOrder = true;
                bool hasFullOrder = (prod.perms.size() == 1);
                for (const Permutation& prodPerm : prod.perms) {
                  if (prodPerm.order == 60) { hasFullOrder = true; }
                  if (prodPerm.order != 60 && prodPerm.order != 30 && prodPerm.order != 20 && prodPerm.order != 12 && prodPerm.order != 1) {
                    hasValidOrder = false;
                    break;
                  }
                }
                if (!hasValidOrder) { continue; }
                if (embeddable && !hasFullOrder) { continue; }

                //Print result
                os << (rep5E && rep3E && rep2zE && rep2xE && hasFullOrder ? " " : "*");
                os << "I  : " << (rep5Ix + 1) << "," << (sym3Ix + 1) << "," << (sym2zIx + 1) << "," << (sym2xIx + 1)
                   << "   " << rep5.id << " x" << sym3.id << " x" << sym2z.id << " x" << sym2x.id << std::endl;
                found = true;
                break;
              }
              if (found) { break; }
            }
            if (found) { break; }
          }
        }
      }
    }
  }
}
