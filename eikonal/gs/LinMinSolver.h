#ifndef _LIN_MIN_SOLVER_H
#define _LIN_MIN_SOLVER_H

#include "Solver.h"

namespace eikonal
{
  class LinMinSolver : public Solver
  {
    /*
      Solver of the Eikonal equation based on linear brent minimization problem.
    */
  public:
    // constructor
    LinMinSolver(const dolfin::FunctionSpace& V);

    // destructor
    virtual ~LinMinSolver();
  
  public:
    static std::string name;

  private:
    // local solver of the Eikonal equation     
    virtual double local_solver(const dolfin::la_index unset_dof,
                            const std::vector<dolfin::la_index>& cell_set_dofs,
                            const dolfin::GenericVector& u_vector,
                            const std::size_t precision);
  };
}

#endif // _LIN_MIN_SOLVER
