#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "gs_connectivity.h"
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>

namespace dolfin
{
  class FunctionSpace;
  class Function;
  class Mesh;
  class GenericVector;
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

    // destructor
    virtual ~Solver();

    // solve the eikonal euqtion |grad(u)| = 1 with values of u fixed at
    // fixed_dofs TODO
    std::size_t solve(dolfin::Function& u,
                      const std::set<dolfin::la_index>& fixed_dofs);
  
  public:
    static std::string name;

  private:
    // initiliaze dof_status TODO, 
    boost::shared_ptr<std::vector<bool> >
    init_dof_status(const dolfin::Function& u,
                    const std::set<dolfin::la_index>& fixed_dofs);

    // initialize unset_dofs, numbers [first,...,last[ exclude fixed TODO
    // fixed_dofs will be sorted inside so no const
    boost::shared_ptr<std::vector<dolfin::la_index> >
    init_unset_dofs(const dolfin::Function& u,
                    const std::set<dolfin::la_index>& fixed_dofs);

    // return dofs in cell_2_dof[cell] that are not dof and whose dof_status
    // is true, i.e set, TODO
    std::vector<dolfin::la_index>
    get_cell_set_dofs(const std::size_t cell, const dolfin::la_index dof) const;

    // get reference points; for now corners of domain assuming square, box
    // domain         TODO
    std::vector<std::vector<double> > 
    get_reference_points(const dolfin::Mesh& mesh) const;

    // local solver of Eikonal equation in the cell
    // compute the approximate value of u_vector in unset_dof from values
    // in set_dofs TODO
    virtual double local_solver(const dolfin::la_index unset_dof,
                             const std::vector<dolfin::la_index>& cell_set_dofs,
                             const dolfin::GenericVector& u_vector) const;

  protected:
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
    boost::shared_ptr<std::vector<bool> > dof_status;

    // indices of all dofs of V, assignment also in solve
    boost::shared_ptr<std::vector<dolfin::la_index> > unset_dofs;

    // first index of dof on the process
    std::size_t offset;
  };
}

#endif // _SOLVER_H_
