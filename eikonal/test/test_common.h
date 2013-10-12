#ifndef _TEST_COMMON_H_
#define _TEST_COMMON_H_

#include <set>

namespace dolfin
{
  class Cell;
  class Mesh;
  template<typename T> class MeshFunction;
}

namespace eikonal
{
  // see it the triangular cell of the mesh is accute
  bool accute_cell(const dolfin::Cell& cell);

  // get all cells that are obtuse by their indices
  std::set<std::size_t> asset_obtuse_cells(const dolfin::Mesh& mesh);

  // get all cells that are obtuse as a mesh function
  dolfin::MeshFunction<bool> obtuse_cells(const dolfin::Mesh& mesh);
}

#endif // _TEST_COMMON_H_
