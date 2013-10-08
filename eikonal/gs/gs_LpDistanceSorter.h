#ifndef _GS_LP_DISTANCE_SORTER_H_
#define _GS_LP_DISTANCE_SORTER_H_

#include <map>
#include <vector>
#include <dolfin/common/types.h>       // la_index
#include "gs_MyIterator.h"

namespace eikonal
{
  class LpDistanceSorter
  {
  /*
    Sort dofs=std::vector<dolfin::la_index> by their l^p distance from
    reference points=std::vector<std::vector<double> >. Dof coordinates are
    looked up in the the map; map[dof] = dof_coordinates.
  */

  public:
    // constructor, set the map
    LpDistanceSorter(const std::map<std::size_t, std::vector<double> >& _map) :
    map(_map) { }

    // sort dofs by their l^p distance from reference points
    void sort(std::vector<dolfin::la_index>& dofs,
              const std::vector<std::vector<double> >& ref_points,
              const std::size_t p);
   
    // get iterator to dofs sorted in (reverse) order from k-th point
    MyIterator<dolfin::la_index> get(const std::size_t k, bool reverse) const;

    // functor for sorting
    bool operator()(const dolfin::la_index& i, const dolfin::la_index& j) const;

  private:
    // set ref point for sorting
    void set_ref_point(const std::vector<double>& _point)
    { ref_point = _point; }

    // set sort l^p norm
    void set_p(const std::size_t _p) { p = _p; }

  private:
    // the map
    const std::map<std::size_t, std::vector<double> >& map;

    // dofs sorted accorded to distance from k-th points
    std::map<std::size_t, std::vector<dolfin::la_index> > sorted_dofs;

    // reference point currently used by functor for sorting
    std::vector<double> ref_point;

    // p of the l^p norm currently used by functor for sorting
    std::size_t p;
  };
}

#endif // _GS_LP_DISTANCE_SORTER_H_
