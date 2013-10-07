#include "Problem.h"
#include "Seeder.h"
#include <set>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Point.h>

using namespace dolfin;

namespace eikonal
{
  Problem::Problem(const eikonal::Seeder& _seeder) : seeder(_seeder) { }
  //---------------------------------------------------------------------------

  void Problem::init(std::vector<dolfin::la_index> fixed_dofs,
                     dolfin::Function& u) const
  {
    // get the points from seeder
    std::vector<Point> points;
    seeder.seed(points, 1000);

    // get the intersected cell
    std::set<std::size_t> cells;
    u.function_space()->mesh()->intersected_cells(points, cells);

    // push to fixed_dofs and modify u
    fixed_dofs.clear();
    boost::shared_ptr<const GenericDofMap> dofmap(u.function_space()->dofmap());
    std::set<std::size_t>::const_iterator cell;
    for(cell = cells.begin(); cell != cells.end(); cell++)
    {
      // fixed dofs
      std::vector<la_index> cell_dofs = dofmap->cell_dofs(*cell);
      fixed_dofs.insert(fixed_dofs.end(), cell_dofs.begin(), cell_dofs.end());

      // get coordinates of dofs and compute distance for setting u
      boost::multi_array<double, 2> dof_coordinates;
      dofmap->tabulate_coordinates(dof_coordinates,
                                   UFCCell(Cell(*u.function_space()->mesh(),
                                           *cell)));

      std::vector<la_index>::const_iterator dof = cell_dofs.begin();
      int index = 0;
      for( ; dof != cell_dofs.end(); dof++, index++)
      {
        std::vector<double> dof_x(dof_coordinates.shape()[1]);// dof_coordinate
        for(std::size_t i = 0; i < dof_x.size(); i++)
        {
          dof_x[i] = dof_coordinates[index][i]; // VERY AWKWARD!!! FIXME
        }
        u.vector()->setitem(*dof, seeder.distance(dof_x));
      }
    }
  }
  //---------------------------------------------------------------------------

  void Problem::exact_solution(dolfin::Function& u) const
  {

  }
  //---------------------------------------------------------------------------
}



