#ifndef _LIN_NEWTON_SOLVER_H
#define _LIN_NEWTON_SOLVER_H

#include "Solver.h"

namespace eikonal
{
  class LinNewtonSolver : public Solver
  {
    /*
      Solver of the Eikonal equation based on linear newton minimization problem.
    */
  public:
    // constructor
    LinNewtonSolver(const dolfin::FunctionSpace& V);

    // destructor
    virtual ~LinNewtonSolver();
  
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

#endif // _LIN_NEWTON_SOLVER
