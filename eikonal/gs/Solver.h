#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "gs_connectivity.h"
#include <vector>

namespace dolfin
{
  class FunctionSpace;
  class Function;
}

namespace eikonal
{
  class Solver
  {
    /*
      Preliminary global solver based on fast sweeping and geometrical solver.
    */
  public:
    // constructor
    Solver(const dolfin::FunctionSpace& V);

    // solve the eikonal euqtion |grad(u)| = 1 if u fixed at fixed_dofs
    void solve(dolfin::Function& u,
               std::vector<dolfin::la_index>& fixed_dofs); 

  public: // protected


  private:
    // return dofs in cell_2_dof[cell] that are not dof and whose dof_status
    // is true, i.e set
    std::vector<dolfin::la_index> get_cell_set_dofs(const std::size_t cell,
                                                const dolfin::la_index dof,
                                const std::vector<bool>& dof_status) const;
 
    // let dof_status be std::map<la_index, bool>, change its values corrensponding
    // to indices in fixed_dofs to true
    //update_dof_status(std::vector<bool>& dof_status,
    //                  const std::vector<dolfin::la_index>& fixed_dofs) const;

  private:
    // initiliaze _dof_status
    


  private:
    // map between cell=size_t and dofs it contains=set([la_index])
    const eikonal::t_smap_la cell_2_dof;

    // inverse of cell_2_dof
    const eikonal::la_vmap_t dof_2_cell;

    // map dof=size_t to its coordinates
    const eikonal::t_vmap_d dof_2_coordinate;

    // reference to function space for possible later use
    const dolfin::FunctionSpace& V;

    // status of (some) dofs of V, this should be assigned in solve(...),
    // when it becomes clear which dofs we need
    std::vector<bool>* dof_status;

    // indices of all dofs of V, assignment also in solve
    std::vector<dolfin::la_index>* unset_dofs;
  };
}

#endif // _SOLVER_H_
