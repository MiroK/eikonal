#include "Problem.h"
#include "Seeder.h"
#include <vector>
#include <set>
#include <cstdio>         //
#include <cstdlib>        // random
#include <ctime>          //
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshFunction.h>
#include "la/la_common.h"

#include <dolfin/plot/plot.h>

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
    

  void Problem::init(std::set<dolfin::la_index>& fixed_dofs,
                     dolfin::Function& u, dolfin::Function& du_dx,
                     dolfin::Function& du_dy) const
  {
    // get the points from seeder
    std::vector<Point> points;
    seeder.seed(points, 1000);

    // get the intersected cell
    std::set<std::size_t> cells;
    u.function_space()->mesh()->intersected_cells(points, cells);

    // push to fixed_dofs and modify u in fixed_dofs, other dofs set as far
    const double far = u.function_space()->mesh()->hmax()*
                       u.function_space()->mesh()->num_cells();
    boost::shared_ptr<GenericVector> u_vector(u.vector());
    *u_vector += far;
    // vectors for the gradient components
    boost::shared_ptr<GenericVector> du_dx_vector(du_dx.vector());
    boost::shared_ptr<GenericVector> du_dy_vector(du_dy.vector());
    
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

        std::vector<double> grad_u;
        seeder.gradient(dof_x, grad_u);
        du_dx_vector->setitem(*dof, grad_u[0]);
        du_dy_vector->setitem(*dof, grad_u[1]);
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

  void Problem::exact_solution(dolfin::Function& u, dolfin::Function& du_dx,
                               dolfin::Function& du_dy) const
  {
    // call to values
    exact_solution(u);

    // get the gradients, //TODO checks for same space for alll function

    // get coordinates of all dofs
    boost::shared_ptr<const GenericDofMap> dofmap(u.function_space()->dofmap());
    std::vector<double> dof_coordinates =
    dofmap->tabulate_all_coordinates(*u.function_space()->mesh());

    std::size_t num_dofs = u.function_space()->dim();
    std::size_t dim = u.function_space()->mesh()->geometry().dim();

    std::vector<double> dx_values(num_dofs);
    std::vector<double> dy_values(num_dofs);
    for(std::size_t i = 0; i < num_dofs; i++)
    {
      const std::vector<double> point(&dof_coordinates[i*dim],
                                      &dof_coordinates[i*dim] + dim);
      std::vector<double> gradient;
      seeder.gradient(point, gradient);
      dx_values[i] = gradient[0];
      dy_values[i] = gradient[1];
    }
    du_dx.vector()->set_local(dx_values);
    du_dy.vector()->set_local(dy_values);
  }
  //---------------------------------------------------------------------------
  
  dolfin::MeshFunction<std::size_t>
  Problem::get_band(dolfin::Function& u, const std::size_t width) const
  {
    // get the points from seeder
    std::vector<Point> points;
    seeder.seed(points, 1000);
    
    // get the intersected cells;
    boost::shared_ptr<const Mesh> mesh(u.function_space()->mesh()); 
    std::set<std::size_t> cells;
    mesh->intersected_cells(points, cells);

    // mark the intersected cells 1
    const std::size_t dim = mesh->topology().dim();
    MeshFunction<std::size_t> band(*mesh, dim);
    band.set_all(0);

    for(std::set<std::size_t>::const_iterator cell = cells.begin();
        cell != cells.end(); cell++)
    {
      band[*cell] = 1;
    }

    // now mark the neighbors, get the cell-cell connectivity
    mesh->init(dim, dim);
    std::set<std::size_t>::const_iterator marked_cell;
    for(std::size_t level = 0; level < width; level++)
    {
      std::set<std::size_t> marked_cells;
      for(CellIterator cell(*mesh); !cell.end(); ++cell)
      {
        if(band[*cell] == 1)
        {
          std::size_t num_next = cell->num_entities(dim);
          const std::size_t* next = cell->entities(dim);
          marked_cells.insert(next, next+num_next);
        }
      }
      marked_cell = marked_cells.begin();
      for( ; marked_cell != marked_cells.end(); marked_cell++)
      {
        band[*marked_cell] = 1;
      }
    }

    return band;
  }
  //---------------------------------------------------------------------------

  void
  Problem::get_ref_points(const std::size_t n_ref_points,
                          std::vector<std::vector<double> >& ref_points) const
  {
    // get the points from seeder
    std::vector<Point> points;
    seeder.seed(points, 1000);

    ref_points.erase(ref_points.begin(), ref_points.end());
    // check how points seeder actually provided
    std::size_t n_points = points.size();
    if(n_points < n_ref_points)
    {
      for(std::size_t i = 0; i < n_points; i++)
      {
        const double* x = points[i].coordinates();
        ref_points.push_back(std::vector<double>(x, x + 2));
      }
    }
    else
    {
      std::set<std::size_t> used_points;
      srand(time(NULL));
      for(std::size_t i = 0; i < n_ref_points; i++)
      {
        std::size_t ref_point;
        do
        {
          ref_point = rand() % n_points;
        } while(used_points.find(ref_point) != used_points.end());
        used_points.insert(ref_point);

        const double* x = points[ref_point].coordinates();
        ref_points.push_back(std::vector<double>(x, x + 2));
      }
    }
  }
  //---------------------------------------------------------------------------

  std::string Problem::name() const { return seeder.name; }
  //---------------------------------------------------------------------------
}
