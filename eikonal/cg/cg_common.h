#ifndef _CG_COMMON_H_
#define _CG_COMMON_H_

/*
  Functions used or computing with polygons and polyhedra.
*/

#include <vector>

namespace eikonal
{
  
  // get the berycenter of poly(gon/hedra) with vertices in R^dim
  std::vector<double> barycenter(const std::vector<double>& vertices,
                                 const std::size_t dim);
}
#endif //_CG_COMMON_H_
