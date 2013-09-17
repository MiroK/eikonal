#include "connectivity.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <sstream>

using namespace dolfin;

//-----------------------------------------------------------------------------
tMapLa cell_to_dof(const dolfin::FunctionSpace& V)
{
  boost::shared_ptr<const Mesh> mesh = V.mesh();
  boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();

  tMapLa _cell_to_dof;
  for(CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    std::size_t cell_index = cell->index();
    _cell_to_dof[cell_index] = dofmap->cell_dofs(cell_index);
  }

  return _cell_to_dof;
}

//-----------------------------------------------------------------------------
laMapT dof_to_cell(const tMapLa& _cell_to_dof)
{
  laMapT _dof_to_cell;
  tMapLa::const_iterator item;
  for(item = _cell_to_dof.begin(); item != _cell_to_dof.end(); item++)
  { 
    std::vector<la_index>::const_iterator dof = item->second.begin();
    for( ; dof != item->second.end(); dof++)
    {
      _dof_to_cell[*dof].push_back(item->first); // push_back(cell)
    }
  }

  return _dof_to_cell;
}

//-----------------------------------------------------------------------------
tMapT cell_to_facet(const dolfin::FunctionSpace& V)
{
  boost::shared_ptr<const Mesh> mesh = V.mesh();

  tMapT _cell_to_facet;
  for(CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    std::size_t cell_index = cell->index();
    for(EdgeIterator facet(*cell); !facet.end(); ++facet)
      _cell_to_facet[cell_index].push_back(facet->index());
  }

  return _cell_to_facet;
}

//-----------------------------------------------------------------------------
tMapT facet_to_cell(const tMapT& _cell_to_facet)
{
  tMapT _facet_to_cell;
  tMapT::const_iterator item;
  for(item = _cell_to_facet.begin(); item != _cell_to_facet.end(); item++)
  {
    std::vector<std::size_t>::const_iterator facet = item->second.begin();
    for(; facet != item->second.end(); facet++)
    {
      _facet_to_cell[*facet].push_back(item->first);
    }
  }
  
  return _facet_to_cell;
}

//-----------------------------------------------------------------------------
tMapT facet_to_dof(const dolfin::FunctionSpace& V)
{
  boost::shared_ptr<const Mesh> mesh = V.mesh();
  boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();

  tMapT _facet_to_dof;
  for(CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    std::size_t i = 0;
    std::vector<la_index> cell_dofs = dofmap->cell_dofs(cell->index());
    for(EdgeIterator facet(*cell); !facet.end(); ++facet)
    {
      std::vector<std::size_t> facet_dofs;
      dofmap->tabulate_facet_dofs(facet_dofs, i);
      std::vector<std::size_t>::iterator dof = facet_dofs.begin();
      for( ; dof != facet_dofs.end(); dof++)
      {
        *dof = cell_dofs[*dof];
      }
      
      _facet_to_dof[facet->index()] = facet_dofs;
      i++;
    }
  }

  return _facet_to_dof;
}

//-----------------------------------------------------------------------------
tMapT dof_to_facet(const tMapT& _facet_to_dof)
{
  tMapT _dof_to_facet;
  tMapT::const_iterator item;
  for(item = _facet_to_dof.begin(); item != _facet_to_dof.end(); item++)
  {
    std::vector<std::size_t>::const_iterator dof = item->second.begin();
    for(; dof != item->second.end(); dof++)
    {
      _dof_to_facet[*dof].push_back(item->first);
    }
  }
  
  return _dof_to_facet;
}
