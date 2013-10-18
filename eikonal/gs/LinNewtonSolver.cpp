#include "LinNewtonSolver.h"
#include "ls/ls_minimize.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericVector.h>

namespace eikonal
{
  std::string LinNewtonSolver::name = std::string("linear_2d_newton");

  LinNewtonSolver::LinNewtonSolver(const dolfin::FunctionSpace& V) : Solver(V) { }
  //----------------------------------------------------------------------------

  LinNewtonSolver::~LinNewtonSolver(){ }
  //----------------------------------------------------------------------------

  double 
  LinNewtonSolver::local_solver(const dolfin::la_index unset_dof,
                             const std::vector<dolfin::la_index>& cell_set_dofs,
                             const dolfin::GenericVector& u_vector)
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
    
    std::size_t n_calls;
    double u_new = linear_newton_2d(u_point, u_value, k_points, k_values,
                                    n_calls);
    
    if(n_calls > max_calls)
    {
      max_calls = n_calls;
    }
    if(n_calls < min_calls)
    {
      min_calls = n_calls; 
    }
    return u_new;
  }
  //----------------------------------------------------------------------------
}
