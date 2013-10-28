#ifndef _HERMITE_SOLVER_H_
#define _HERMITE_SOLVER_H_

#include "Solver.h"
#include <utility>

namespace eikonal
{
  class Sorter;

  class HermiteSolver : public Solver
  {
    /* Solver based on fast sweeping method with 
    Hermite polynomial reconstruction of the front*/

  public:
    // constructor, V is CG1
    HermiteSolver(const dolfin::FunctionSpace& V);

    // destructor
    virtual ~HermiteSolver();

    // p controls the norm used in sorter
    // ref points can be empty then the reference points are computed
    std::size_t solve(dolfin::Function& u,
                      dolfin::Function& du_dx,
                      dolfin::Function& du_dy,
                      const std::set<dolfin::la_index>& fixed_dofs,
                      const std::size_t precision,
                      const std::size_t p,
                      std::vector<std::vector<double> >& ref_points);
 
    // unset dofs are set by Sorter
    std::size_t solve(dolfin::Function& u,
                      dolfin::Function& du_dx,
                      dolfin::Function& du_dy,
                      const std::set<dolfin::la_index>& fixed_dofs,
                      const std::size_t precision,
                      const Sorter& sorter); 
  public:
    static std::string name;
  

  private:
    // local solver of the eikonal equation, with dofs ordered by distance
    // from 4 reference - edges of the mesh
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
