#include "FMMSolver.h"
#include "ls/ls_geometric.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <algorithm>

#include "la/la_common.h"

using namespace dolfin;

namespace eikonal
{
  FMMSolver::FMMSolver(const dolfin::FunctionSpace& _V) :
                                  cell_2_dof(cell_to_dof(_V)),
                                  dof_2_cell(dof_to_cell(cell_2_dof)),
                                  dof_2_coordinate(dof_to_coordinate(_V)),
                                  dof_2_dof(dof_to_dof(cell_2_dof, dof_2_cell)),
                                  offset(0)
  {
    // TODO assertions on V, ideally before all the mapping build
  }
  //---------------------------------------------------------------------------

  FMMSolver::~FMMSolver(){ }
  //---------------------------------------------------------------------------
  
  void FMMSolver:: solve(dolfin::Function& u, std::set<dolfin::la_index>& fixed)
  {
    // initialize the vector of close dofs, unigue dofs connected to fixed
    // that are not fixed or close themselves
    
    std::cout << "fixed_dofs: "; print(fixed);
    
    std::vector<la_index> close;
    update_close(fixed, fixed, close);
    
    // FMM
    boost::shared_ptr<GenericVector> u_vector = u.vector();
    Compare compare(*u_vector);
    std::vector<la_index> update = close;
    while(not close.empty())
    {
      std::cout << "close: "; print(close);
      // update distance for dofs in update
      std::vector<la_index>::const_iterator dof = update.begin();
      for( ; dof != update.end(); dof++)
      {
        double u_old = (*u_vector)[*dof];
        std::vector<std::size_t> cells = dof_2_cell.find(*dof)->second;

        std::cout << "\tdof: " << *dof << std::endl;
        std::cout << "\tcells: "; print(cells);

        std::vector<std::size_t>::const_iterator cell = cells.begin();
        for( ; cell != cells.end(); cell++)
        {
          std::vector<la_index>
          cell_fixed_dofs = get_cell_fixed_dofs(*cell, *dof, fixed);
          std::cout << "\t\tdofs: " << * cell; print(cell_fixed_dofs);

          double u_ = this->local_solver(*dof, cell_fixed_dofs, *u_vector);
          std::cout << "\t\tcell produced: " << u_ << std::endl;
          if(u_ < u_old)
          {
            u_old = u_;
            std::cout << "\t\t\tfound minimizer\n";
          }
        }
        u_vector->setitem(*dof, u_old);
      }

      // get the dof with the smallest distance and remove it from closest
      std::make_heap(close.begin(), close.end(), compare);
      std::pop_heap(close.begin(), close.end(), compare);
      la_index closest = close.back();
      std::cout << "closest is " << closest << std::endl;
      close.pop_back();

      // add closest to fixed and his neighbors to close
      fixed.insert(closest);
      std::cout << "fixed now "; print(fixed);

      // add the neighbors of closest tha are not fixed and close to close
      // mark all the neighbors for next round of distance computation
      update = std::vector<la_index>(dof_2_dof.find(closest)->second.begin(),
                                     dof_2_dof.find(closest)->second.end());
      std::cout << "dofs to update "; print(update);

      std::set<la_index> dofs;
      dofs.insert(closest);
      update_close(dofs, fixed, close);
    }
  }
  //---------------------------------------------------------------------------

  std::vector<dolfin::la_index>
  FMMSolver::get_cell_fixed_dofs(const std::size_t cell, 
                           const dolfin::la_index dof,
                           const std::set<dolfin::la_index>& fixed_dofs) const
  {
    std::vector<la_index> cell_set_dofs;
    std::set<la_index>::const_iterator cell_dof =
    cell_2_dof.find(cell)->second.begin();
    for( ; cell_dof != cell_2_dof.find(cell)->second.end(); cell_dof++)
    {
      if(*cell_dof != dof and fixed_dofs.find(*cell_dof) != fixed_dofs.end())
      {
        cell_set_dofs.push_back(*cell_dof);
      }
    }

    return cell_set_dofs;
  }
  //---------------------------------------------------------------------------

  void FMMSolver::update_close(const std::set<dolfin::la_index>& dofs,
                               const std::set<dolfin::la_index> fixed,
                               std::vector<dolfin::la_index>& close)
  {
    std::set<la_index>::const_iterator _dof = dofs.begin();
    for(; _dof != dofs.end(); _dof++)
    {
      std::set<la_index> connected = dof_2_dof.find(*_dof)->second;
      std::set<la_index>::const_iterator dof = connected.begin();
      for(; dof != connected.end(); dof++)
      {
        const bool in_fixed = 
        std::find(fixed.begin(), fixed.end(), *dof) != fixed.end();

        const bool in_close = 
        std::find(close.begin(), close.end(), *dof) != close.end();
        
        if(not in_fixed and not in_close)
        {
          close.push_back(*dof);
        }
      }
    }
  }
  //---------------------------------------------------------------------------


  double FMMSolver::local_solver(const dolfin::la_index close_dof,
                          const std::vector<dolfin::la_index>& cell_fixed_dofs,
                          const dolfin::GenericVector& u_vector)
  {
    const std::size_t size = cell_fixed_dofs.size();
    double new_value;
    if(size == 2 or size == 1)
    {
      // prepare data for the solver with ls interface
      const std::vector<double>
      u_point = dof_2_coordinate.find(close_dof)->second;
      const double u_value = u_vector[close_dof];
      
      const std::size_t g_dim = 2;
      double _k_values[size];
      double _k_points[g_dim*size];
      for(std::size_t fixed_dof = 0; fixed_dof < size; fixed_dof++)
      {
        // set value of known dofs
        _k_values[fixed_dof] = u_vector[cell_fixed_dofs[fixed_dof]];

        // set coodinates of known dofs
        std::vector<double> dof_coordinate =
        dof_2_coordinate.find(cell_fixed_dofs[fixed_dof])->second;
        _k_points[g_dim*(fixed_dof)+0] = dof_coordinate[0];
        _k_points[g_dim*(fixed_dof)+1] = dof_coordinate[1];
      }
      const std::vector<double> k_points(_k_points, _k_points + g_dim*size);
      const std::vector<double> k_values(_k_values, _k_values + 2);
     
      if(size == 2)
      {
        new_value = linear_geometric_2d(u_point, u_value, k_points, k_values);
        std::cout << "\t\tlocal_solver 2 value " << new_value;
      }
      else
      {
        new_value = linear_extrapolate(u_point, k_points, k_values[0]); 
      }
    }
    else
    {
      new_value = u_vector[close_dof];
    }

    return new_value;
  }
  //---------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------

namespace eikonal
{
  Compare::Compare(const dolfin::GenericVector& _value) : value(_value) { }
  //---------------------------------------------------------------------------
    
  bool Compare::operator()
  (const dolfin::la_index& i, const dolfin::la_index& j) const
  {
    return value[i] > value[j];
  }
  //---------------------------------------------------------------------------
}
