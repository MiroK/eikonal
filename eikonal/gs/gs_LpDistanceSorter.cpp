#include "gs_LpDistanceSorter.h"
#include "la/la_loop.h"               // - of vectors and norm
#include <algorithm>
#include <cassert>

using namespace dolfin;

namespace eikonal
{
  void LpDistanceSorter::sort(std::vector<dolfin::la_index>& dofs,
                              const std::vector<std::vector<double> >& ref_points,
                              const std::size_t p)
  {
    set_p(p);
    for(std::size_t k = 0; k < ref_points.size(); k++)
    {
      set_ref_point(ref_points[k]);
      std::sort(dofs.begin(), dofs.end(), *this);
      sorted_dofs[k] = dofs;
    }
  }
  //-----------------------------------------------------------------------------
   
  MyIterator<dolfin::la_index>
  LpDistanceSorter::get(const std::size_t k, bool reverse) const
  {
    assert(k < sorted_dofs.size());
    return MyIterator<la_index>(sorted_dofs.find(k)->second, reverse);
  }
  //-----------------------------------------------------------------------------

  bool LpDistanceSorter::
  operator()(const dolfin::la_index& i, const dolfin::la_index& j) const
  {
    return norm(map.find(i)->second - ref_point, p) <
           norm(map.find(j)->second - ref_point, p);
  }
  //---------------------------------------------------------------------------
}
