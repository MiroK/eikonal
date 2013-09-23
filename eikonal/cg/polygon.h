#ifndef _POLYGON_H_
#define _POLYGON_H_

/*
  Functions for computations with polygons.
*/

#include <vector>

// create regular n-gon inscribed into circle with center and radius
// center has size 2 and vertices
std::vector<double>
polygon_generate(const std::size_t n, const std::vector<double>& center,
                 const double radius);

// given vertices of polygon (flattened [v0_x, v0_y, v1_x, v1_y, ... ])
// compute vertices that are on boundary of polygon
// t > 1 specifies number of vertices per edge of polygon
std::vector<double>
polygon_boundary(const std::vector<double>& vertices, const std::size_t t = 100);

// compute area of polygon specified by vertices
double polygon_area(const std::vector<double>& vertices);

// compute sum of area of triangles formed by point and pair of vertices
double polygon_area_wrt(const std::vector<double>& vertices,
                        const std::vector<double>& point);

// check if point is inside the polygon
bool polygon_inside(const std::vector<double>& vertices,
                    const std::vector<double>& point);

#endif // _POLYGON_H_
