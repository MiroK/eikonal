#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

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

// types of map that are build
typedef std::map<std::size_t, std::vector<dolfin::la_index> > tMapLa; 
typedef std::map<dolfin::la_index, std::vector<std::size_t> > laMapT;
typedef std::map<std::size_t, std::vector<std::size_t> > tMapT;

//--------------------MAP BUILDING---------------------------------------------

// build a map: map[cell] = [dofs in cell]; size_t -> vector(la_index)
tMapLa cell_to_dof(const dolfin::FunctionSpace& V);

// build a map: map[dof] = [cells with dof]; la_index -> vector(size_t)
laMapT dof_to_cell(const tMapLa& _cell_to_dof);

// build a map: map[cell] = [facets of the cell]; size_t -> vector(size_t)
tMapT cell_to_facet(const dolfin::FunctionSpace& V);

// build a map: map[facet] = [cells that share facet]; size_t -> vector(size_t)
tMapT facet_to_cell(const tMapT& _cell_to_facet);

// build a map: map[facet] = [dofs on the facet]; size_t -> vector(size_t)
tMapT facet_to_dof(const dolfin::FunctionSpace& V);

// build a map: map[dof] = [facets that share dof]; size_t -> vector(size_t)
tMapT dof_to_facet(const tMapT& _facet_to_dof);

//--------------------CONVENIENCE----------------------------------------------

template<typename T, typename U>
std::string print(const T& map, const std::string map_name);

template<typename I, typename O, typename U>
O invert_map(const I& input_map);

//-------------------implementations-------------------------------------------
template<typename T, typename U>
std::string print(const T& map, const std::string map_name)
{
  std::ostringstream out;
  out << "\t" << map_name << std::endl;
  
  typename T::const_iterator item;
  for(item = map.begin(); item != map.end(); item++)
  {
    std::size_t key = item->first;
    out << key << " : ";
    typename std::vector<U> values = item->second;
    typename std::vector<U>::const_iterator value;
    for(value = values.begin(); value != values.end(); value++)
    {
      out << *value << " ";
    }
    out << std::endl;
  }

  return out.str();
}


template<typename I, typename O, typename U>
O invert_map(const I& input_map)
{
  O inverse_map;
  typename I::const_iterator item;
  for(item = input_map.begin(); item != input_map.end(); item++)
  { 
    typename std::vector<U>::const_iterator value = item->second.begin();
    for( ; value != item->second.end(); value++)
    {
      inverse_map[*value].push_back(item->first); // push_back(cell)
    }
  }

  return inverse_map;
}
#endif // _CONNECTIVITY_H_
