#pragma once
#include "types.h"
#include <vector>

enum class SymType {
  None,
  Mirror,
  Chiral,
  Spiegel,
  Dihedral,
  Tetrahedral,
  Octahedral,
  Icosahedral,
};

class Permutation {
public:
  Permutation() : order(0) {}
  Permutation(const std::vector<int>& _cycle) : order(0), cycle(_cycle) {}

  //Representative of the group.
  int rep() const { return cycle[0]; }
  //Index of the or element -1 if not found.
  int indexOf(int element) const;

  //Make all cycles the same length and save original order
  bool Normalize(int max_order);

  //Apply the symmetry to the actual 3d points or planes
  bool Apply(VectorXf& x, SymType type) const;

  //Number of permutations before repeating.
  int order;
  std::vector<int> cycle;
};
