#ifndef _CG_POLYGON_H_
#define _CG_POLYGON_H_

/*
  Functions for computations with polygons. 
*/

#include <vector>

namespace eikonal
{
  // create regular n-gon inscribed into circle with center and radius
  // returns flattened vector containing vertices [v0_x, v0_y, v1_x, v1_y, ... ]
  // center and vertices are 2D
  std::vector<double> g_vertices(const std::size_t n,
                                 const std::vector<double>& center,
                                 const double radius);

  // given vertices of polygon (flattened )
  // compute poins that are on boundary of polygon
  // pps > 1 specifies number of vertices per edge of polygon
  std::vector<double> boundary_points(const std::vector<double>& vertices,
                                      const std::size_t pps = 100);

  // compute area of polygon specified by vertices
  double area(const std::vector<double>& vertices);

  // compute sum of areas of triangles formed by point and pair of vertices
  // i.e. sum of barycentric coordinates
  double area_wrt(const std::vector<double>& vertices,
                  const std::vector<double>& point);

  // check if point is inside the polygon
  bool g_inside(const std::vector<double>& vertices,
                const std::vector<double>& point);

  // get the barycenter of polygon specified by vertices
  std::vector<double> g_barycenter(const std::vector<double>& vertices);

  // see if triangle with vertices [v0x, v0y, v1x...] is accute,
  // has angle greater than pi/2
  bool accute_triangle(std::vector<double>& vertices);
}

#endif // _CG_POLYGON_H_
