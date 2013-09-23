#include "polygon.h"
#include "polyhedron.h"
#include <cmath>
#include <stdexcept>

using namespace std;

std::vector<double>
polygon_generate(const std::size_t n, const std::vector<double>& center,
                 const double radius)
{
  if(center.size() == 2)
  {
    std::vector<double> vertices(2*n);
    const double c_x = center[0];
    const double c_y = center[1];
    const double d_theta = 2.*M_PI/n;
    for(size_t i = 0; i < n; i++)
    {
      vertices[2*i] = c_x + radius*cos(i*d_theta);
      vertices[2*i + 1] = c_y + radius*sin(i*d_theta);
    }
    return vertices;
  }
  else
    throw length_error("polygon.cpp polygon_generate(...) center size must be 2");
}

//-----------------------------------------------------------------------------

std::vector<double>
polygon_boundary(const std::vector<double>& vertices, const std::size_t t)
{
}

//-----------------------------------------------------------------------------

double polygon_area(const std::vector<double>& vertices)
{
  const size_t n_vertices = vertices.size()/2;
  const std::vector<double> center = barycenter(vertices, n_vertices);
  return polygon_area_wrt(vertices, center);
}

//-----------------------------------------------------------------------------

double polygon_area_wrt(const std::vector<double>& vertices,
                        const std::vector<double>& point)
{
  double area = 0;
}

//-----------------------------------------------------------------------------

bool polygon_inside(const std::vector<double>& vertices,
                    const std::vector<double>& point)
{

}
