#ifndef _CG_POLYHEDRON_H_
#define _CG_POLYHEDRON_H_

/*
  Functions for computations with with tetrahedrons.
*/

#include <vector>

namespace eikonal
{

  // compute barycenter of polyhedron specified by vertices
  // 4 vertices of size 3 flattened 
  std::vector<double> t_barycenter(const std::vector<double>& vertices);

  // volume of tetrahedon given by vertices
  double volume(const std::vector<double>& vertices);

  // volume with respect to point, i.e. sum of barycentric coordinates
  double volume_wrt(const std::vector<double>& vertices,
                    const std::vector<double>& point);

  // see if point is inside tetrahedron
  double t_inside(const std::vector<double>& vertices,
                  const std::vector<double>& point);

  // return ppe points per every edge of the tetrahedon given by vertices
  std::vector<double> edge_points(const std::vector<double>& vertices,
                                  const std::size_t ppe=100);

}

#endif // _CG_POLYHEDRON_H_
