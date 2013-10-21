#include "HermiteSolver.h"
#include "gs/gs_LpDistanceSorter.h"
#include "ls/ls_minimize.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>

using namespace dolfin;

namespace eikonal
{
  std::string HermiteSolver::name = std::string("hermite_2d");
  
  HermiteSolver::HermiteSolver(const dolfin::FunctionSpace& V) : Solver(V) { }
  //---------------------------------------------------------------------------
  
  HermiteSolver::~HermiteSolver(){ }
  //---------------------------------------------------------------------------

  std::size_t HermiteSolver::solve(dolfin::Function& u,
                                   dolfin::Function& du_dx,
                                   dolfin::Function& du_dy,
                                   const std::set<dolfin::la_index>& fixed_dofs,
                                   const std::size_t precision)
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
    boost::shared_ptr<GenericVector> du_dx_vector = du_dx.vector();
    boost::shared_ptr<GenericVector> du_dy_vector = du_dy.vector();

    Function v(V);                // solution from previous weep
    boost::shared_ptr<GenericVector> v_vector = v.vector();
    *v_vector = *u_vector;

    // sweep
    std::size_t k = 0; // counter of unique sweeps;
    while(true)
    {
      k++;
      std::size_t ref_point = ((k-1)%n_uniq_sweeps)/2;
      bool reverse_flag = (k-1)%2 ?  true : false;
      // apply local solver to unset_dofs in order given by ref_point
     
      for(MyIterator<la_index> unset_dof = sorter.get(ref_point, reverse_flag);
          !unset_dof.end(); ++unset_dof)
      {
        double u_old = (*u_vector)[*unset_dof];
        double du_dx_old = (*du_dx_vector)[*unset_dof];
        double du_dy_old = (*du_dy_vector)[*unset_dof];
        std::vector<std::size_t> cells = dof_2_cell.find(*unset_dof)->second;

        std::vector<std::size_t>::const_iterator cell = cells.begin();
        for( ; cell != cells.end(); cell++)
        {
          std::vector<la_index>
          cell_set_dofs = get_cell_set_dofs(*cell, *unset_dof);

          std::pair<double, std::vector<double> > u_gradu = 
          this->local_solver(*unset_dof, cell_set_dofs, *u_vector, 
                             *du_dx_vector, *du_dy_vector, precision);
          double u_ = u_gradu.first;
          if(u_ < u_old)
          {
            u_old = u_;
            (*dof_status)[*unset_dof - offset] = true;
            // also set the gradient
            du_dx_old = (u_gradu.second)[0];
            du_dy_old = (u_gradu.second)[1];
          }
        }
        u_vector->setitem(*unset_dof, u_old);
        du_dx_vector->setitem(*unset_dof, du_dx_old);
        du_dy_vector->setitem(*unset_dof, du_dy_old);
      }

      // sweep over; check convergence as |u-v| in L infty norm
      *v_vector -= *u_vector;
      v_vector->abs();
      if(v_vector->max() < DOLFIN_EPS)
      {
        return k;
      }
      else
      {
        *v_vector = *u_vector;
      }
    }
  }
  //---------------------------------------------------------------------------
  
  std::pair<double, std::vector<double> > 
  HermiteSolver::local_solver(const dolfin::la_index unset_dof,
                              const std::vector<dolfin::la_index>& cell_set_dofs,
                              const dolfin::GenericVector& u_vector,
                              const dolfin::GenericVector& du_dx_vector,
                              const dolfin::GenericVector& du_dy_vector,
                              const std::size_t precision)
  {
    std::vector<double> grad_u_C(2);
    
    if(cell_set_dofs.size() != 2)
    {
      return std::make_pair<double, std::vector<double> >(u_vector[unset_dof],
                                                                    grad_u_C);
    }

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
    
    std::size_t n_calls;
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
    return std::make_pair<double, std::vector<double> >(t_ft.second, grad_u_C);
  }
  //--------------------------------------------------------------------------
}
