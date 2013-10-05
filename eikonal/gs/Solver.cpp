#include "Solver.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>

using namespace dolfin;

namespace eikonal
{
  Solver::Solver(const dolfin::FunctionSpace& _V) :
                                      cell_2_dof(cell_to_dof(_V)),
                                      dof_2_cell(dof_to_cell(cell_2_dof)),
                                      dof_2_coordinate(dof_to_coordinate(_V)),
                                      V(_V)
  {
    //TODO check assertions on function space
    dof_status = NULL;
    unset_dofs = NULL;
  }
  //---------------------------------------------------------------------------

  void Solver::solve(dolfin::Function& u,
                     std::vector<dolfin::la_index>& fixed_dofs)
  {
    // change dof_status
    /*std::size_t dim  = u.function_space()->dim();
    std::vector<bool> _dof_status(dim, false);
    for(std::size_t i = 0; i < fixed_dofs.size(); i++)
      _dof_status[fixed_dofs[i]] = true;
    
    dof_status = &_dof_status;

    // change dofs_to_set
    std::vector<la_index> _dofs_to_set(dim);
    for(std::size_t i = 0; i < dim; i++)
    {
      _dofs_to_set[i] = i;
    }

    std::sort(fixed_dofs.begin(), fixed_dofs.end());
    std::vector<la_index>::const_reverse_iterator fixed_it;
    std::vector<la_index>::iterator 
    for(dof_it = fixed_dofs.rbegin(); dof_it != fixed_dofs.rend(); dof_it++)
    {
      _dofs_to_set.erase(dof_it); 
    }

    dofs_to_set = &_dofs_to_set;

    // create reference points and order dofs_to_set

    // sweep
    std::size_t iter_count = 0;
    boost::shared_ptr<GenericVector> u_vector = u.vector();
    while(true)
    {
      iter_count++;
      std::size_t k = (iter_count-1) % 8; // we assume for reference points
                         // with each point there are two sweeps associated
      for()


    }*/
  }
  //---------------------------------------------------------------------------

  std::vector<dolfin::la_index>
  Solver::get_cell_set_dofs(std::size_t cell, dolfin::la_index dof,
                            const std::vector<bool>& dof_status) const
  {
    std::vector<la_index> cell_set_dofs;
    std::set<la_index>::const_iterator cell_dof =
    cell_2_dof.find(cell)->second.begin();
    for( ; cell_dof != cell_2_dof.find(cell)->second.end(); cell_dof++)
    {
      std::cout << *cell_dof << std::endl;
      if(dof != *cell_dof and dof_status[*cell_dof])
      {
        cell_set_dofs.push_back(*cell_dof);
      }
    }

    if(cell_set_dofs.size() != 2)
    {
      cell_set_dofs.clear();
    }
    return cell_set_dofs;
  }
  //---------------------------------------------------------------------------
}
