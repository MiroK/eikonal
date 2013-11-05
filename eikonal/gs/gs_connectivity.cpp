#include "gs_connectivity.h"
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <sstream>

using namespace dolfin;

namespace eikonal
{
  t_smap_la cell_to_dof(const dolfin::FunctionSpace& V)
  {
    boost::shared_ptr<const Mesh> mesh = V.mesh();
    boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();

    t_smap_la _cell_to_dof;
    for(CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      std::size_t cell_index = cell->index();
      std::vector<la_index> dofs = dofmap->cell_dofs(cell_index);
      _cell_to_dof[cell_index] = std::set<la_index>(dofs.begin(), dofs.end());
    }

    return _cell_to_dof;
  }
  //-----------------------------------------------------------------------------

  la_vmap_t dof_to_cell(const t_smap_la& _cell_to_dof)
  {
    return invert_map<t_smap_la, la_vmap_t, la_index>(_cell_to_dof);
  }
  //---------------------------------------------------------------------------

  la_smap_la dof_to_dof(const t_smap_la& _cell_to_dof,
                        const la_vmap_t& _dof_to_cell)
  {
    la_smap_la _dof_to_dof;
    la_vmap_t::const_iterator dc_entry = _dof_to_cell.begin();
    for( ; dc_entry != _dof_to_cell.end(); dc_entry++)
    {
      std::size_t dof = dc_entry->first;
      
      // insert dofs of cells with dof;
      std::vector<la_index>
      cells(dc_entry->second.begin(), dc_entry->second.end());
      std::vector<la_index>::const_iterator cell = cells.begin();
      for( ; cell != cells.end(); cell++)
      {
        _dof_to_dof[dof].insert(_cell_to_dof.at(*cell).begin(),
                                _cell_to_dof.at(*cell).end());
      }

      // remove the dof
      _dof_to_dof[dof].erase(dof);
    }

    return _dof_to_dof;
  }
  //-----------------------------------------------------------------------------
  
  t_smap_t cell_to_facet(const dolfin::FunctionSpace& V)
  {
    boost::shared_ptr<const Mesh> mesh = V.mesh();

    t_smap_t _cell_to_facet;
    for(CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      std::size_t cell_index = cell->index();
      for(FacetIterator facet(*cell); !facet.end(); ++facet)
        _cell_to_facet[cell_index].insert(facet->index());
    }

    return _cell_to_facet;
  }

  //-----------------------------------------------------------------------------
  t_vmap_t facet_to_cell(const t_smap_t& _cell_to_facet)
  {
    return invert_map<t_smap_t, t_vmap_t, std::size_t >(_cell_to_facet);
  }

  //-----------------------------------------------------------------------------
  t_smap_t facet_to_dof(const dolfin::FunctionSpace& V)
  {
    boost::shared_ptr<const Mesh> mesh = V.mesh();
    boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();

    t_smap_t _facet_to_dof;
    for(CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      std::size_t i = 0; // 3 edges for triangle, 4 faces for tetrahedra
      std::vector<la_index> cell_dofs = dofmap->cell_dofs(cell->index());
      // facets could be obtained from cell_to_facet but this way we have
      // consistent interaface
      for(FacetIterator facet(*cell); !facet.end(); ++facet)
      {
        std::vector<std::size_t> facet_dofs;
        dofmap->tabulate_facet_dofs(facet_dofs, i);
        std::vector<std::size_t>::iterator dof = facet_dofs.begin();
        for( ; dof != facet_dofs.end(); dof++)
        {
          *dof = cell_dofs[*dof];
        }
        
        _facet_to_dof[facet->index()] = std::set<std::size_t>(facet_dofs.begin(),
                                                              facet_dofs.end());
        i++;
      }
    }

    return _facet_to_dof;
  }

  //-----------------------------------------------------------------------------
  t_vmap_t dof_to_facet(const t_smap_t& _facet_to_dof)
  {
    return invert_map<t_smap_t, t_vmap_t, std::size_t >(_facet_to_dof);
  }

  //-----------------------------------------------------------------------------
  t_vmap_d dof_to_coordinate(const dolfin::FunctionSpace& V)
  {
    boost::shared_ptr<const Mesh> mesh = V.mesh();
    boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();
    std::vector<double> all_coordinates = dofmap->tabulate_all_coordinates(*mesh);

    std::size_t gdim = mesh->geometry().dim();
    std::size_t n_dofs = all_coordinates.size()/gdim;
    
    t_vmap_d _dof_to_coordinate;
    std::vector<double>::const_iterator start = all_coordinates.begin();
    for(std::size_t i = 0; i < n_dofs; i++)
    {
      _dof_to_coordinate[i].assign(start + i*gdim, start + (i+1)*gdim);
    }

    return _dof_to_coordinate;
  }
}
