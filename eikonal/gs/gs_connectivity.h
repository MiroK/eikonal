#ifndef _GS_CONNECTIVITY_H_
#define _GS_CONNECTIVITY_H_

/*
  Functions for obtaining cell-dof-cell, cell-facet-cell, facet-dof-facet
  connectivity.
*/

#include <dolfin/common/types.h>
#include <boost/shared_ptr.hpp>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <sstream>

namespace dolfin // Forward declarations
{
  class FunctionSpace;
}

namespace eikonal
{
  // types of map that are build
  typedef std::map<std::size_t, std::set<dolfin::la_index> > t_smap_la; 
  typedef std::map<dolfin::la_index, std::vector<std::size_t> > la_vmap_t;
  typedef std::map<std::size_t, std::set<std::size_t> > t_smap_t; 
  typedef std::map<std::size_t, std::vector<std::size_t> > t_vmap_t;
  typedef std::map<std::size_t, std::vector<double> > t_vmap_d;

  //--------------------MAP BUILDING---------------------------------------------

  // build a map: map[cell] = [dofs in cell]; size_t -> set(la_index)
  t_smap_la cell_to_dof(const dolfin::FunctionSpace& V);

  // build a map: map[dof] = [cells with dof]; la_index -> vector(size_t)
  la_vmap_t dof_to_cell(const t_smap_la& _cell_to_dof);

  // build a map: map[cell] = [facets of the cell]; size_t -> set(size_t)
  t_smap_t cell_to_facet(const dolfin::FunctionSpace& V);

  // build a map: map[facet] = [cells that share facet]; size_t -> vector(size_t)
  t_vmap_t facet_to_cell(const t_smap_t& _cell_to_facet);

  // build a map: map[facet] = [dofs on the facet]; size_t -> set(size_t)
  t_smap_t facet_to_dof(const dolfin::FunctionSpace& V);

  // build a map: map[dof] = [facets that share dof]; size_t -> vector(size_t)
  t_vmap_t dof_to_facet(const t_smap_t& _facet_to_dof);

  // build a map: map[dof] = [dof coordinates]
  t_vmap_d dof_to_coordinate(const dolfin::FunctionSpace& V);

  // print content of map of type T with values of type U
  template<typename T, typename U>
  std::string print_map(const T& map, const std::string map_name);

  // invert the map of type I with values of type std::set<U> into map of type O
  template<typename I, typename O, typename U>
  O invert_map(const I& map);

  //-------------------implementations-------------------------------------------
  template<typename T, typename U>
  std::string print_map(const T& map, const std::string map_name)
  {
    std::ostringstream out;
    out << "\t" << map_name << std::endl;
    
    typename T::const_iterator item;
    for(item = map.begin(); item != map.end(); item++)
    {
      std::size_t key = item->first;
      out << key << " : ";
      U values = item->second;
      typename U::const_iterator value;
      for(value = values.begin(); value != values.end(); value++)
      {
        out << *value << " ";
      }
      out << std::endl;
    }

    return out.str();
  }

  template<typename I, typename O, typename U>
  O invert_map(const I& map)
  {
    O inverse_map;
    typename I::const_iterator item;
    for(item = map.begin(); item != map.end(); item++)
    { 
      typename std::set<U>::const_iterator value = item->second.begin();
      for( ; value != item->second.end(); value++)
      {
        inverse_map[*value].push_back(item->first);
      }
    }

    return inverse_map;
  }
}

#endif // _GS_CONNECTIVITY_H_
