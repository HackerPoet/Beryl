#pragma once
#include "types.h"
#include "permutation.h"
#include <iostream>
#include <vector>
#include <unordered_set>

class Symmetry {
public:
  Symmetry() { Reset(); }

  //Clears the symmetry to default values.
  void Reset();

  //Make all permutations the same length (with repetition) and save orders.
  //Must be called after symmetry is fully constructed.
  bool Normalize();

  //Set the type of symmetry for optimization.
  void SetType(SymType _type) { type = _type; }

  //Is symmetry uninitialized
  bool Empty() const { return perms.empty(); }

  //Permute an element using the symmetry
  int Permute(int element, int shift=1) const;
  std::vector<int> Permute(const std::vector<int>& nums, int shift=1) const;
  static std::vector<int> IdentityPerm(size_t size);
  static bool IsIdentityPerm(const std::vector<int>& nums);

  //Load symmetry from automorphism file.
  static bool Load(Symmetry& sym, Symmetry& dsym, const Faces& tris, const char* fname, int line_num, bool zero_indexed=false);
  static bool Load(Symmetry& sym, Symmetry& dsym, const Faces& tris, const char* fname, const std::vector<int>& line_nums, bool zero_indexed=false);

  //Multiply symmetries to form a new one.
  static Symmetry Multiply(const std::vector<Symmetry>& syms);

  //Find the dual symmetry
  Symmetry Dual(const Faces& tris) const;

  //Apply the symmetry to the actual 3d points or planes
  bool Apply(VectorXf& x) const;

  //Create a unique id for this permutation based on how far it sends each vertex
  std::string MakeId(const Faces& tris, const std::map<Edge, int>& dist_map);
  static std::map<Edge, int> MakeDistMap(const Faces& tris);

  //Test if a set of vertices has this symmetry and return the name and axis
  std::string Test(const Verts3D& x) const;

  //Get symmetry enum from string
  static SymType GetType(char c);
  static int GetProdSize(SymType symType);

  SymType type;
  int order;
  int numElements;
  std::string id;
  std::unordered_set<int> orders;
  std::vector<Permutation> perms;
};

class SymmetryList {
public:
  void Reset();

  bool Load(const Faces& tris, const char* fname, bool zero_indexed = false);

  void Analyze(const Faces& tris, bool embeddable, std::ostream& os = std::cout) const;

  std::vector<Symmetry> syms;
  std::vector<Symmetry> symDuals;
  std::map<int, std::vector<int>> symsByOrder;
  std::map<std::string, std::vector<int>> symsById;
};
