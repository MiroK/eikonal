#ifndef _HERMITE_SOLVER_H_
#define _HERMITE_SOLVER_H_

#include "Solver.h"
#include <utility>

namespace eikonal
{
  class HermiteSolver : public Solver
  {
    /* Solver based on fast sweeping method with 
    Hermite polynomial reconstruction of the front*/

  public:
    // constructor, V is CG1
    HermiteSolver(const dolfin::FunctionSpace& V);

    // destructor
    virtual ~HermiteSolver();

    std::size_t solve(dolfin::Function& u,
                      dolfin::Function& du_dx,
                      dolfin::Function& du_dy,
                      const std::set<dolfin::la_index>& fixed_dofs,
                      const std::size_t precision);
  
  public:
    static std::string name;
  

  private:
    // local solver of the eikonal equation
    std::pair<double, std::vector<double> >
    local_solver(const dolfin::la_index unset_dof,
                 const std::vector<dolfin::la_index>& cell_set_dofs,
                 const dolfin::GenericVector& u_vector,
                 const dolfin::GenericVector& du_dx_vector,
                 const dolfin::GenericVector& du_dy_vector,
                 const std::size_t precision);
  };
}

#endif // _HERMITE_SOLVER_H
