#ifndef _FMM_HERMITE_H_
#define _FMM_HERMITE_H_

#include "FMMSolver.h"
#include <utility>

namespace eikonal
{
  class Sorter;

  class FMMHermite : public FMMSolver
  {
    /*
     * Fast marching method based on the hermite interpolation
     */

    public:
      // constructor
      FMMHermite(const dolfin::FunctionSpace& V) : FMMSolver(V) { }

      // destructor
      virtual ~FMMHermite() { }

      // solve the eikonal equation
      std::size_t solve(dolfin::Function& u,
                        dolfin::Function& du_dx,
                        dolfin::Function& du_dy,
                        std::set<dolfin::la_index>& fixed,
                        const std::size_t precision);
      
      // these are dummies for compatibility
      std::size_t solve(dolfin::Function& u,
                        dolfin::Function& du_dx,
                        dolfin::Function& du_dy,
                        std::set<dolfin::la_index>& fixed_dofs,
                        const std::size_t precision,
                        const std::size_t p,
                        std::vector<std::vector<double> >& ref_points,
                        const std::size_t max_sweep=100);
   
      // these are dummies for compatibility
      std::size_t solve(dolfin::Function& u,
                        dolfin::Function& du_dx,
                        dolfin::Function& du_dy,
                        std::set<dolfin::la_index>& fixed_dofs,
                        const std::size_t precision,
                        const Sorter& sorter,
                        const std::size_t max_sweep=100); 
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

#endif // _FMM_HERMITE_H_
