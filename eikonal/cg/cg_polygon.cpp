#include "cg_polygon.h"
#include "cg_common.h"     // barycenter
#include "la/la_loop.h"    // vector operations
#include "la/la_common.h"  // compare
#include <cmath>        
#include <cassert>        
#include <algorithm>

using namespace std;

namespace eikonal
{
  std::vector<double> g_vertices(const std::size_t n,
                                 const std::vector<double>& center,
                                 const double radius)
  {
    assert(center.size() == 2);
    vector<double> vertices(2*n);
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
  //---------------------------------------------------------------------------

  std::vector<double> boundary_points(const std::vector<double>& vertices,
                                      const std::size_t pps)
  {
    assert(pps > 1);
    
    const size_t dim = 2;
    assert(vertices.size()%dim == 0);

    const size_t n_vertices = vertices.size()/dim;
    
    const double dt = 1./(pps-1);
    const size_t n_points = n_vertices*(pps - 1);
    vector<double> boundary_points;
    boundary_points.reserve(dim*n_points);
    vector<double>::iterator boundary_point = boundary_points.begin();
    
    for(size_t i = 0;  i < n_vertices; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[((i+1)%n_vertices)*dim];
      const vector<double> B(_B, _B + dim);

      for(size_t j = 0; j < (pps - 1); j++)
      {
        const double t = j*dt;
        const vector<double> point = A + t*(B-A);
        boundary_points.insert(boundary_point, point.begin(), point.end());
        boundary_point += dim;
      }
    }
    return boundary_points;
  }
  //---------------------------------------------------------------------------
  
  double area(const std::vector<double>& vertices)
  {
    // assert handled in area_wrt
    const vector<double> x0 = g_barycenter(vertices);
    return area_wrt(vertices, x0);
  }
  //---------------------------------------------------------------------------

  double area_wrt(const std::vector<double>& vertices,
                  const std::vector<double>& point)
  {
    assert(point.size() == 2);
    
    const size_t dim = 2;
    assert(vertices.size()%dim == 0);
    
    const size_t n_vertices = vertices.size()/dim;
    double _area=0;
    for(size_t i = 0;  i < n_vertices; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[((i+1)%n_vertices)*dim];
      const vector<double> B(_B, _B + dim);

      const vector<double> x = eikonal::cross(B-point, A-point);
      _area += eikonal::norm(x, 2);
    }
  
    return _area/2.;
  }
  //---------------------------------------------------------------------------

  bool g_inside(const std::vector<double>& vertices,
                const std::vector<double>& point)
  {
    // assert handled in area_wrt
    return eikonal::compare(area(vertices), area_wrt(vertices, point));
  }
  //---------------------------------------------------------------------------


  std::vector<double> g_barycenter(const std::vector<double>& vertices)
  {
    assert(vertices.size()%2 == 0);
    return eikonal::barycenter(vertices, 2);
  }
  //---------------------------------------------------------------------------

  bool accute_triangle(std::vector<double>& vertices)
  {
    assert(vertices.size() == 6);

    // get edge lengths
    std::vector<double> edge_sizes(3);
    for(std::size_t i = 0; i < 3; i++)
    {
      std::vector<double> P(&vertices[2*i], &vertices[2*i]+2);
      std::vector<double> Q(&vertices[2*((i+1)%3)], &vertices[2*((i+1)%3)]+2);
      edge_sizes[i] = norm(P-Q, 2);
    }

    // pick the longest edge, largest angle is opposite to it
    std::size_t i = std::distance(edge_sizes.begin(),
                        std::max_element(edge_sizes.begin(), edge_sizes.end()));
    std::vector<double> C(&vertices[2*i], &vertices[2*i]+2);
    std::vector<double> A(&vertices[2*((i+1)%3)], &vertices[2*((i+1)%3)]+2);
    std::vector<double> B(&vertices[2*((i+2)%3)], &vertices[2*((i+2)%3)]+2);
    
    // cos of largest angle
    const double cosine = dot(A-C, B-C)/edge_sizes[(i+1)%3]/edge_sizes[(i+2)%3];
    return not (cosine > 1);
  }
  //----------------------------------------------------------------------------
}
