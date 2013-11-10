#ifndef _FMM_SOLVER_H_
#define _FMM_SOLVER_H_

#include "gs_connectivity.h"
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>

namespace dolfin
{
  class Function;
  class GenericVector;
}

namespace eikonal
{
  class FMMSolver
  {
    /*
      Preliminary version of global EIkonal solver based on Fast marching
      method. As a local solver, geometrical solver is used.
    */
  public:
    // constructor
    FMMSolver(const dolfin::FunctionSpace& V);

    // destructor
    virtual ~FMMSolver();

    // solve the eikonal euqtion |grad(u)| = 1 with values of u fixed at
    // fixed_dofs, precision is the parameter used with iteratative local solver
    // convergence tolerance = std::num_limits<double>::digits/precision
    std::size_t solve(dolfin::Function& u,
                std::set<dolfin::la_index>& fixed,
                std::size_t precision,
                std::size_t compat_dummy=0);

  protected:
    // using cell->dof mapping extract from cell all the dofs that are not dof
    // and are fixed
    std::vector<dolfin::la_index>
    get_cell_fixed_dofs(const std::size_t cell,
                        const dolfin::la_index dof,
                        const std::set<dolfin::la_index>& fixed_dofs) const;
    
    // using dof->dof mapping exctract all dofs that are that are connected
    // to dofs and are not fixed and close already TODO
    void update_close(const std::set<dolfin::la_index>& dofs,
                      const std::set<dolfin::la_index> fixed,
                      std::vector<dolfin::la_index>& close);
    
  private:
    // apply the local solver
    double
    local_solver(const dolfin::la_index close_dof,
                 const std::vector<dolfin::la_index>& cell_fixed_dofs,
                 const dolfin::GenericVector& u_vector);

  protected:
    // map between cell=size_t and dofs it contains=set([la_index])
    const eikonal::t_smap_la cell_2_dof;

    // inverse of cell_2_dof
    const eikonal::la_vmap_t dof_2_cell;

    // map dof=size_t to its coordinates
    const eikonal::t_vmap_d dof_2_coordinate;

    // connectivity of dofs
    const eikonal::la_smap_la dof_2_dof;

    // parallel offset
    std::size_t offset;
    
  public:
    static std::string name;
    std::size_t min_calls;
    std::size_t max_calls;

  };
  
  //---------------------------------------------------------------------------
  
  class Compare
  {
    /*
     * Comparison operator for the min_heap
     */
  public:
    // constructor
    Compare(const dolfin::GenericVector& _value);
    
    // comparison
    bool operator()(const dolfin::la_index& i, const dolfin::la_index& j) const;

  private:  
    // vector of value used for sort
    const dolfin::GenericVector& value;

  };
}

#endif // _FMM_SOLVER_H_
