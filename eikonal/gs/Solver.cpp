#include "Solver.h"
#include "ls/ls_geometric.h"
#include "gs/gs_LpDistanceSorter.h"
#include "la/la_loop.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/common/constants.h>
#include <algorithm>
#include <cmath>

#include "la/la_common.h"

using namespace dolfin;

namespace eikonal
{
  std::string Solver::name = std::string("linear_2d_geometric");

  Solver::Solver(const dolfin::FunctionSpace& _V) :
                                      cell_2_dof(cell_to_dof(_V)),
                                      dof_2_cell(dof_to_cell(cell_2_dof)),
                                      dof_2_coordinate(dof_to_coordinate(_V)),
                                      V(_V),
                                      offset(0),
                                      min_calls(1E6),
                                      max_calls(0)
  {
    // TODO assertions on V, ideally before all the mapping build
  }
  //---------------------------------------------------------------------------

  Solver::~Solver(){ }
  //---------------------------------------------------------------------------

  boost::shared_ptr<std::vector<bool> >
  Solver::init_dof_status(const dolfin::Function& u,
                          const std::set<dolfin::la_index>& fixed_dofs)
  {
    //[status of dof_first, ..., dof_second)
    std::pair<std::size_t, std::size_t> u_range = u.vector()->local_range();

    if(offset == -1) // set offset if it wasn't set before
    {
      offset = u_range.first;
    }
    
    // initialize all to false for start
    boost::shared_ptr<std::vector<bool> >
    _dof_status(new std::vector<bool>(u_range.second - offset, false));

    // look through fixed_dofs, setting corresponding entries in dof_status to
    // true, make sure fixed_dofs is on the process
    std::set<la_index>::const_iterator fixed_dof = fixed_dofs.begin();
    for( ; fixed_dof != fixed_dofs.end(); fixed_dof++)
    {
      if(*fixed_dof >= offset and *fixed_dof < u_range.second)
      {
        (*_dof_status)[*fixed_dof - offset] = true;
      }
    }

    return _dof_status;
  }
  //---------------------------------------------------------------------------

  boost::shared_ptr<std::vector<dolfin::la_index> >
  Solver::init_unset_dofs(const dolfin::Function& u,
                          const std::set<dolfin::la_index>& fixed_dofs)
  {
    std::pair<std::size_t, std::size_t> u_range = u.vector()->local_range();

    if(offset == -1) // set offset if it wasn't set before
    {
      offset = u_range.first;
    }

    // put all dofs into unset at first
    boost::shared_ptr<std::vector<la_index> >
    _unset_dofs(new std::vector<la_index>(u_range.second - offset));
    for(std::size_t i = 0; i < _unset_dofs->size(); i++)
    {
      (*_unset_dofs)[i] = offset + i;
    }

    // remove the fixed dofs
    std::set<la_index>::const_reverse_iterator
    fixed_dof = fixed_dofs.rbegin();
    for( ; fixed_dof != fixed_dofs.rend(); fixed_dof++)
    {
      if(*fixed_dof >= offset and *fixed_dof < u_range.second)
      {
        _unset_dofs->erase(_unset_dofs->begin() + (*fixed_dof - offset));
      }
    }

    return _unset_dofs;
  }
  //---------------------------------------------------------------------------
  
  std::vector<dolfin::la_index>
  Solver::get_cell_set_dofs(std::size_t cell, dolfin::la_index dof) const
  {
    std::vector<la_index> cell_set_dofs;
    std::set<la_index>::const_iterator cell_dof =
    cell_2_dof.find(cell)->second.begin();
    for( ; cell_dof != cell_2_dof.find(cell)->second.end(); cell_dof++)
    {
      if(dof != *cell_dof and (*dof_status)[*cell_dof])
      {
        cell_set_dofs.push_back(*cell_dof);
      }
    }

    return cell_set_dofs;
  }
  //---------------------------------------------------------------------------


  std::vector<std::vector<double> > 
  Solver::get_reference_points(const dolfin::Mesh& mesh) const
  {
    std::size_t dim = mesh.geometry().dim();
    std::size_t num_vertices = mesh.num_vertices();
    std::vector<double> coordinates = mesh.coordinates();

    // search for minima and maxima of in every direction
    std::vector<double> min_max(2*dim);
    for(std::size_t d = 0; d < dim; d++)
    {
      double min = coordinates[d], max = coordinates[d];
      for(std::size_t i = 1; i < num_vertices; i++)
      {
        double current = coordinates[i*dim + d];
        
        if(current > max) max = current;

        if(current < min) min = current;

      }
      min_max[2*d] = min;
      min_max[2*d + 1] = max;
    }

    // create ref. points
    std::vector<std::vector<double> > ref_points((std::size_t)pow(2, dim));
    for(std::size_t d = 0; d < dim; d++)
    {
      double d_min = min_max[2*d];
      double d_max = min_max[2*d + 1];
      for(std::size_t i = 0; i < ref_points.size(); i++)
      {
        double e = (i/((std::size_t)pow(2, dim - (d+1))))%2 ? d_max : d_min;
        ref_points[i].push_back(e);
      }
    }
    return ref_points;
  }
  //---------------------------------------------------------------------------
  
  double
  Solver::local_solver(const dolfin::la_index unset_dof,
                       const std::vector<dolfin::la_index>& cell_set_dofs,
                       const dolfin::GenericVector& u_vector,
                       const std::size_t precision)
  {
    if(cell_set_dofs.size() != 2)
    {
      return u_vector[unset_dof];
    }

    // prepare data for the solver with ls interface
    const std::vector<double>
    u_point = dof_2_coordinate.find(unset_dof)->second;
    const double u_value = u_vector[unset_dof];
    
    double _k_points[4];
    for(std::size_t set_dof = 0; set_dof < 2; set_dof++)
    {
      std::vector<double>
      dof_coordinate = dof_2_coordinate.find(cell_set_dofs[set_dof])->second;
      _k_points[2*(set_dof)+0] = dof_coordinate[0];
      _k_points[2*(set_dof)+1] = dof_coordinate[1];
    }
    const std::vector<double> k_points(_k_points, _k_points + 4);
    double _k_values[2] = {u_vector[cell_set_dofs[0]],
                           u_vector[cell_set_dofs[1]]};
    const std::vector<double> k_values(_k_values, _k_values + 2);
    
    // geoemtric does not use precision!
    double new_value = linear_geometric_2d(u_point, u_value, k_points, k_values);
    return new_value;
  }
  //---------------------------------------------------------------------------
  
  
  std::size_t Solver::solve(dolfin::Function& u,
                            const std::set<dolfin::la_index>& fixed_dofs,
                            const std::size_t precision,
                            const std::size_t max_sweep)
  {
    // initialization
    
    dof_status = init_dof_status(u, fixed_dofs);
    
    unset_dofs = init_unset_dofs(u, fixed_dofs);
    
    std::vector<std::vector<double> > ref_points = get_reference_points(*V.mesh());
   
    LpDistanceSorter sorter(dof_2_coordinate);
    // use L^2 norm to sort unset_dofs by their distance from ref_points;
    sorter.sort(*unset_dofs, ref_points, 2); 

    std::size_t n_points = ref_points.size();
    std::size_t n_uniq_sweeps = n_points*2;
    
    boost::shared_ptr<GenericVector> u_vector = u.vector();
    Function v(V);                // solution from previous weep
    boost::shared_ptr<GenericVector> v_vector = v.vector();
    *v_vector = *u_vector;

    // sweep
    std::size_t k = 0; // counter of unique sweeps;
    bool converged = false;
    while(not converged and k < max_sweep)
    {
      k++;
      std::size_t ref_point = ((k-1)%n_uniq_sweeps)/2;
      bool reverse_flag = (k-1)%2 ?  true : false;
      // apply local solver to unset_dofs in order given by ref_point
     
      for(MyIterator<la_index> unset_dof = sorter.get(ref_point, reverse_flag);
          !unset_dof.end(); ++unset_dof)
      {
        double u_old = (*u_vector)[*unset_dof];
        std::vector<std::size_t> cells = dof_2_cell.find(*unset_dof)->second;

        std::vector<std::size_t>::const_iterator cell = cells.begin();
        for( ; cell != cells.end(); cell++)
        {
          std::vector<la_index>
          cell_set_dofs = get_cell_set_dofs(*cell, *unset_dof);

          double u_ = this->local_solver(*unset_dof, cell_set_dofs, *u_vector,
                                         precision);
          if(u_ < u_old)
          {
            u_old = u_;
            (*dof_status)[*unset_dof - offset] = true;
          }
        }
        u_vector->setitem(*unset_dof, u_old);
      }

      // sweep over; check convergence as |u-v| in L infty norm
      *v_vector -= *u_vector;
      v_vector->abs();
      if(v_vector->max() < DOLFIN_EPS)
      {
        converged = true;
      }
      else
      {
        *v_vector = *u_vector;
      }
    }
    return k;
  }
  //---------------------------------------------------------------------------
}
