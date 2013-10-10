#include "Problem.h"
#include "Seeder.h"
#include <vector>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/MeshGeometry.h>

#include "la/la_common.h"

using namespace dolfin;

namespace eikonal
{
  Problem::Problem(const eikonal::Seeder& _seeder) : seeder(_seeder) { }
  //---------------------------------------------------------------------------

  void Problem::init(std::set<dolfin::la_index>& fixed_dofs,
                     dolfin::Function& u) const
  {
    // get the points from seeder
    std::vector<Point> points;
    seeder.seed(points, 1000);

    // get the intersected cell
    std::set<std::size_t> cells;
    u.function_space()->mesh()->intersected_cells(points, cells);

    std::cout << "intersected cells :"; print(cells);

    // push to fixed_dofs and modify u in fixed_dofs, other dofs set as far
    const double far = u.function_space()->mesh()->hmax()*
                       u.function_space()->mesh()->num_cells();
    boost::shared_ptr<GenericVector> u_vector(u.vector());
    *u_vector += far;

    fixed_dofs.clear();
    boost::shared_ptr<const GenericDofMap> dofmap(u.function_space()->dofmap());
    std::set<std::size_t>::const_iterator cell;
    for(cell = cells.begin(); cell != cells.end(); cell++)
    {
      // fixed dofs
      std::vector<la_index> cell_dofs = dofmap->cell_dofs(*cell);
      fixed_dofs.insert(cell_dofs.begin(), cell_dofs.end());

      // get coordinates of dofs and compute distance for setting u
      boost::multi_array<double, 2> dof_coordinates;
      dofmap->tabulate_coordinates(dof_coordinates,
                                   UFCCell(Cell(*u.function_space()->mesh(),
                                           *cell)));
      // len of position vector of dof
      std::size_t dim = dof_coordinates.shape()[1]; 

      std::vector<la_index>::const_iterator dof = cell_dofs.begin();
      int index = 0;
      for( ; dof != cell_dofs.end(); dof++, index++)
      {
        std::vector<double> dof_x(&dof_coordinates[index][0],
                                  &dof_coordinates[index][0] + dim);
        
        u_vector->setitem(*dof, seeder.distance(dof_x));
      }
    }
  }
  //---------------------------------------------------------------------------

  void Problem::exact_solution(dolfin::Function& u) const
  {
    // get coordinates of all dofs
    boost::shared_ptr<const GenericDofMap> dofmap(u.function_space()->dofmap());
    std::vector<double> dof_coordinates =
    dofmap->tabulate_all_coordinates(*u.function_space()->mesh());

    std::size_t num_dofs = u.function_space()->dim();
    std::size_t dim = u.function_space()->mesh()->geometry().dim();

    std::vector<double> values(num_dofs);
    for(std::size_t i = 0; i < num_dofs; i++)
    {
      const std::vector<double> point(&dof_coordinates[i*dim],
                                      &dof_coordinates[i*dim] + dim);
      values[i] = seeder.distance(point);
    }
    u.vector()->set_local(values);
  }
  //---------------------------------------------------------------------------

  std::string Problem::name() const { return seeder.name; }
  //---------------------------------------------------------------------------
}
