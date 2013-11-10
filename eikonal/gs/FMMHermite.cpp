#include "FMMHermite.h"
#include "ls/ls_geometric.h"
#include "ls/ls_minimize.h"
#include "la/la_loop.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <algorithm>

#include "la/la_common.h"

using namespace dolfin;

namespace eikonal
{
  std::string FMMHermite::name =  std::string("fmm_hermite");

  std::size_t FMMHermite::solve(dolfin::Function& u,
                                dolfin::Function& du_dx,
                                dolfin::Function& du_dy,
                                std::set<dolfin::la_index>& fixed,
                                const std::size_t precision)
  {
    // initialize the vector of close dofs, unigue dofs connected to fixed
    // that are not fixed or close themselves
    
    std::cout << "fixed_dofs: "; print(fixed);
    
    std::vector<la_index> close;
    update_close(fixed, fixed, close);
    
    // FMM
    boost::shared_ptr<GenericVector> u_vector = u.vector();
    boost::shared_ptr<GenericVector> du_dx_vector = du_dx.vector();
    boost::shared_ptr<GenericVector> du_dy_vector = du_dy.vector();
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
        double du_dx_old = (*du_dx_vector)[*dof];
        double du_dy_old = (*du_dy_vector)[*dof];

        std::vector<std::size_t> cells = dof_2_cell.find(*dof)->second;

        std::cout << "\tdof: " << *dof << std::endl;
        std::cout << "\tcells: "; print(cells);

        std::vector<std::size_t>::const_iterator cell = cells.begin();
        for( ; cell != cells.end(); cell++)
        {
          std::vector<la_index>
          cell_fixed_dofs = get_cell_fixed_dofs(*cell, *dof, fixed);
          std::cout << "\t\tdofs: " << * cell; print(cell_fixed_dofs);

          std::pair<double, std::vector<double> > u_gradu = 
          this->local_solver(*dof, cell_fixed_dofs, *u_vector, 
                             *du_dx_vector, *du_dy_vector, precision);
          double u_ = u_gradu.first;
          std::cout << "\t\tcell produced: " << u_ << std::endl;
          if(u_ < u_old)
          {
            u_old = u_;
            std::cout << "\t\t\tfound minimizer\n";

            // also set the gradient
            du_dx_old = (u_gradu.second)[0];
            du_dy_old = (u_gradu.second)[1];
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
    return 0;
  }
  
  std::size_t FMMHermite::solve(dolfin::Function& u,
                                dolfin::Function& du_dx,
                                dolfin::Function& du_dy,
                                std::set<dolfin::la_index>& fixed_dofs,
                                const std::size_t precision,
                                const std::size_t p,
                                std::vector<std::vector<double> >& ref_points,
                                const std::size_t max_sweep)
  {
    return solve(u, du_dx, du_dy, fixed_dofs, precision);
  }
  //---------------------------------------------------------------------------
                    
  std::size_t FMMHermite::solve(dolfin::Function& u,
                                dolfin::Function& du_dx,
                                dolfin::Function& du_dy,
                                std::set<dolfin::la_index>& fixed_dofs,
                                const std::size_t precision,
                                const Sorter& sorter,
                                const std::size_t max_sweep)
  {
    return solve(u, du_dx, du_dy, fixed_dofs, precision);
  }
  //---------------------------------------------------------------------------

  std::pair<double, std::vector<double> >
  FMMHermite::local_solver(const dolfin::la_index unset_dof,
                           const std::vector<dolfin::la_index>& cell_set_dofs,
                           const dolfin::GenericVector& u_vector,
                           const dolfin::GenericVector& du_dx_vector,
                           const dolfin::GenericVector& du_dy_vector,
                           const std::size_t precision)
  {
    std::vector<double> grad_u_C(2);
    double new_value;
    const std::size_t size = cell_set_dofs.size();
    
    if(size == 2 or size == 1)
    {
      if(size == 1)
      {
        // extrapolation
        std::vector<double> A = dof_2_coordinate.find(cell_set_dofs[0])->second;
        std::vector<double> C = dof_2_coordinate.find(unset_dof)->second;
        grad_u_C = (C-A)/norm(C-A, 2);

        new_value = linear_extrapolate(C, A, u_vector[cell_set_dofs[0]]); 
      }
      else
      {
        // prepare data for the solver with ls interface
        std::size_t a = cell_set_dofs.at(0);
        std::size_t b = cell_set_dofs.at(1);

        const std::vector<double> C = dof_2_coordinate.find(unset_dof)->second;
        const std::vector<double> A = dof_2_coordinate.find(a)->second;
        const std::vector<double> B = dof_2_coordinate.find(b)->second;
        
        const double u_C = u_vector[unset_dof];
        const double u_A = u_vector[a];
        const double u_B = u_vector[b];

        double _grad_u_A[2];
        _grad_u_A[0] = du_dx_vector[a];
        _grad_u_A[1] = du_dy_vector[a];
        const std::vector<double> grad_u_A(_grad_u_A, _grad_u_A + 2);
        
        double _grad_u_B[2];
        _grad_u_B[0] = du_dx_vector[b];
        _grad_u_B[1] = du_dy_vector[b];
        const std::vector<double> grad_u_B(_grad_u_B, _grad_u_B + 2);
        
        std::size_t n_calls = 0;
        
        std::pair<double, double> t_ft = 
        hermite_newton_2d(A, B, C, u_A, u_B, u_C, grad_u_A, grad_u_B, grad_u_C,
                          n_calls, precision);
       
        if(n_calls > max_calls)
        {
          max_calls = n_calls;
        }
        if(n_calls < min_calls)
        {
          min_calls = n_calls; 
        }

        const double u_max = std::max(u_A, u_B);
        new_value = t_ft.second;
        std::cout << "\t\thermite has " << new_value << std::endl;
        std::cout << "\t\tu_A, u_B " << u_A << " " << u_B << std::endl;
        // zero the gradient and return the old value if results no good
        if(new_value < u_max)
        {
          new_value = u_C;
          grad_u_C.assign(2, 0); // zero the gradient
        }
      }
    }
    else
    {
      // grad is 0 and should not be used outside!
      new_value = u_vector[unset_dof];
    }

    return
    std::make_pair<double, std::vector<double> >(new_value, grad_u_C);
  }
  //---------------------------------------------------------------------------
}
