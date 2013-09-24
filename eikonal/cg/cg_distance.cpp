#include "cg_distance.h"
#include "la/la_loop.h"
#include "la/la_common.h"
#include <algorithm>

using namespace std;

namespace eikonal
{
  double point_point(const std::vector<double>& P,
                     const std::vector<double>& Q)
  {
    return norm(P - Q, 2);
  }
  //---------------------------------------------------------------------------
  
  double point_edge(const std::vector<double>& P,
                       const std::vector<double>& A,
                       const std::vector<double>& B)
  {
    // compute line parameter t of the intersect
    pair<double, double> d_t = point_line(P, A, B, "d");
    const double t = d_t.second;
    if(t < 0 or t > 1) // intersect outside of edge
    {
      return -1;  // distance should be always positeve, so this is a flag
    }
    else
    {
      return d_t.first;
    }
  }
  //---------------------------------------------------------------------------
  
  pair<double, double> 
  point_line(const std::vector<double>& P, const std::vector<double>& A,
                    const std::vector<double>& B, const std::string type)
  {
    assert(type == "sd" or type == "d");
    assert(P.size() == 2 and A.size() == 2 and B.size() == 2);

    // compute line parameter t of the intersect and the intersect
    const double AB = norm(B - A, 2);

    assert(not compare(AB, 0.0)); // |AB-0.0| < LA_EPS is false 

    const double t = dot(P - A, B - A)/AB/AB;
    const vector<double> I = A + t*(B - A);
    if(type == "d") 
    {
      return pair<double, double>(norm(P - I, 2), t);
    }
    else
    {
      // check the sign
      int sign = cross(I-A, P - I)[2] >= 0 ? 1 : -1;
      return pair<double, double>(sign*norm(P - I, 2), t);
    }
  }
  //---------------------------------------------------------------------------
  
  double point_polygon(const std::vector<double>& P,
                       const std::vector<double>& vertices,
                       const std::string type)
  {
    // minimal distance from the vertices and edges of the polygon
    const size_t dim = 2;
    assert(P.size() == 2);
    assert(vertices.size()%dim == 0);
    const size_t n_edges = vertices.size()/dim; // also number of vertices 

    vector<double> distances; // will hold 2*n_edges values
    distances.reserve(2*n_edges);

    //push back distances from verties and edges;
    for(size_t i = 0; i < n_edges; i++)
    {
      const double* _V = &vertices[i*dim];
      const vector<double> V(_V, _V + dim);
      distances.push_back(point_point(P, V));
     
      const double* _W = &vertices[((i+1)%n_edges)*dim];
      const vector<double> W(_W, _W + dim);
      const double d = point_edge(P, V, W);
      if(d != -1)
      {
        distances.push_back(d);
      }
    }
    return *min_element(distances.begin(), distances.end());
  }
  //---------------------------------------------------------------------------
  
  double point_circle(const std::vector<double>& P,
                      const std::vector<double>& center, const double radius,
                      const std::string type)
  {
    assert(P.size() == 2 and center.size() == 2);
    assert(type == "sd" or type == "d");
    if(type == "sd")
    {
      return norm(P - center, 2) - radius;
    }
    else
    {
      return abs(norm(P - center, 2) - radius);
    }
  }
  //---------------------------------------------------------------------------
  
  double point_plane(const std::vector<double>& P,
                     const std::vector<double>& A,
                     const std::vector<double>& B,
                     const std::vector<double>& C,
                     const std::string type)
  {

  }
  //---------------------------------------------------------------------------
  
  double point_3face(const std::vector<double>& P,
                     const std::vector<double>& A,
                     const std::vector<double>& B,
                     const std::vector<double>& C,
                     const std::string type)
  {

  }
  //---------------------------------------------------------------------------
    
  double point_tet(const std::vector<double>& P,
                   const std::vector<double>& vertices,
                   const std::string type)
  {

  }
  //---------------------------------------------------------------------------
  
  double point_sphere(const std::vector<double>& P,
                      const std::vector<double>& center, const double radius,
                      const std::string type)
  {
    assert(P.size() == 3 and center.size() == 3);
    assert(type == "sd" or type == "d");
    if(type == "sd")
    {
      return norm(P - center, 2) - radius;
    }
    else
    {
      return abs(norm(P - center, 2) - radius);
    }
  }
}
