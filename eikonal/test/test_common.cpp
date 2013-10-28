#include "test_common.h"
#include "cg/cg_polygon.h"
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <cassert>

#include <la/la_common.h>

using namespace dolfin;

namespace eikonal
{

  bool accute_cell(const dolfin::Cell& cell)
  {
    assert(cell.dim() == 2);

    const unsigned int* indices = cell.entities(0);
    std::vector<double> vertices(6);
    for(std::size_t i = 0; i < 3; i++)
    {
      vertices[2*i] = cell.mesh().coordinates()[2*indices[i]];
      vertices[2*i+1] = cell.mesh().coordinates()[2*indices[i]+1];
    }

    return accute_triangle(vertices);
  }
  //----------------------------------------------------------------------------

  std::set<std::size_t> asset_obtuse_cells(const dolfin::Mesh& mesh)
  {
    assert(mesh.topology().dim() == 2);
    
    std::set<std::size_t> obtuse_cells;
    for(CellIterator cell(mesh); !cell.end(); ++cell)
    {
      if(not accute_cell(*cell))
      {
        obtuse_cells.insert(cell->index());
      }
    }
    return obtuse_cells;
  }
  //----------------------------------------------------------------------------

  dolfin::MeshFunction<bool> obtuse_cells(const dolfin::Mesh& mesh)
  {
    assert(mesh.topology().dim() == 2);
    MeshFunction<bool> obtuse_cells(mesh, 2);
    obtuse_cells.set_all(false);
    
    for(CellIterator cell(mesh); !cell.end(); ++cell)
    {
      std::cout << cell->index() << " " << accute_cell(*cell) << std::endl;
      if(not accute_cell(*cell))
      {
        obtuse_cells.set_value(cell->index(), true);
      }
    }

    return obtuse_cells;
  }
  //----------------------------------------------------------------------------
}
