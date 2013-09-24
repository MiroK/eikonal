#include "cg_tetrahedron.h"
#include "cg_common.h"     // barycenter
#include "la/la_loop.h"    // vector operations
#include "la/la_common.h"  // compare
#include <cassert>        

using namespace std;

namespace eikonal
{

  std::vector<double> t_barycenter(const std::vector<double>& vertices)
  {
    const size_t dim = 3;
    assert(vertices.size()%dim == 0);
    const size_t n_vertices = vertices.size()/dim;
    assert(n_vertices == 4);
    return eikonal::barycenter(vertices, dim);
  }
  //---------------------------------------------------------------------------

  double volume(const std::vector<double>& vertices)
  {
    return volume_wrt(vertices, t_barycenter(vertices));
  }
  //---------------------------------------------------------------------------

  double volume_wrt(const std::vector<double>& vertices,
                    const std::vector<double>& point)
  {
    assert(point.size() == 3);

    const size_t dim = 3;
    const size_t n_vertices = vertices.size()/dim;
    assert(n_vertices == 4);
    
    double _volume = 0;
    for(size_t i = 0; i < n_vertices; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[((i+1)%n_vertices)*dim];
      const vector<double> B(_B, _B + dim);

      const double* _C = &vertices[((i+2)%n_vertices)*dim];
      const vector<double> C(_C, _C + dim);
  
      _volume += abs(dot(point-A, cross(C-A, B-A)));
    }

    return _volume/6.;
  }
  //---------------------------------------------------------------------------

  double t_inside(const std::vector<double>& vertices,
                  const std::vector<double>& point)
  {
    return eikonal::compare(volume(vertices), volume_wrt(vertices, point));
  }
  //---------------------------------------------------------------------------

  std::vector<double> edge_points(const std::vector<double>& vertices,
                                  const std::size_t ppe)
  {
    assert(ppe > 1); 

    const size_t dim = 3;
    const size_t n_vertices = vertices.size()/dim;
    assert(n_vertices == 4);

    const double dt = 1./(ppe - 1);
    const size_t n_points = n_vertices + (ppe-2)*6;
    vector<double> edge_points;
    edge_points.reserve(dim*n_points);
    vector<double>::iterator edge_point = edge_points.begin();
 
    // first cover vertices by sliding filter: 01, 12, 23, 30
    for(size_t i = 0;  i < n_vertices; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[((i+1)%n_vertices)*dim];
      const vector<double> B(_B, _B + dim);

      for(size_t j = 0; j < (ppe - 1); j++)
      {
        const double t = j*dt;
        const vector<double> point = A + t*(B-A);
        edge_points.insert(edge_point, point.begin(), point.end());
        edge_point += dim;
      }
    }

    // cover remaining edges: 02 13
    for(size_t i = 0;  i < 2; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[(i+2)*dim];
      const vector<double> B(_B, _B + dim);

      for(size_t j = 1; j < (ppe - 1); j++)
      {
        const double t = j*dt;
        const vector<double> point = A + t*(B-A);
        edge_points.insert(edge_point, point.begin(), point.end());
        edge_point += dim;
      }
    }
    return edge_points;
  }
  //---------------------------------------------------------------------------
}


