#include "cg_distance.h"
#include "cg_polygon.h"              // g_inside
#include "cg_tetrahedron.h"          // t_inside
#include "la/la_loop.h"              // vector operations
#include "la/la_common.h"            // compare
#include <algorithm>                 // min_element

using namespace std;

namespace eikonal
{
  double point_point(const std::vector<double>& P,
                     const std::vector<double>& Q)
  {
    return norm(P - Q, 2);
  }
  //---------------------------------------------------------------------------
  
  void point_point_gradient(const std::vector<double>& P,
                            const std::vector<double>& Q,
                            std::vector<double>& gradient)
  {
    const double distance = point_point(P, Q);
    const std::size_t dim = P.size(); 
    
    gradient.resize(dim);
    if(distance > LA_EPS)
    {
      gradient = (P-Q)/distance;
    }
    else
    {
      gradient.assign(dim, 0.0); 
    }
  }
  //---------------------------------------------------------------------------
  
  double point_edge(const std::vector<double>& P,
                       const std::vector<double>& A,
                       const std::vector<double>& B)
  {
    // compute line parameter t of the intersect
    const pair<double, bool> d_inside = point_line(P, A, B, "d");
    const bool inside = d_inside.second;
    if(not inside) // intersect outside of edge
    {
      return -1;  // distance should be always positeve, so this is a flag
    }
    else
    {
      return d_inside.first;
    }
  }
  //---------------------------------------------------------------------------
  
  void point_edge_gradient(const std::vector<double>& P,
                           const std::vector<double>& A,
                           const std::vector<double>& B,
                           std::vector<double>& gradient)
  {
    // get the distance to see if point intersects edge
    double d = point_edge(P, A, B);
    if(d < 0)// no intersect, zero the gradient
    {
      gradient.assign(P.size(), 0.);
    }
    else
    {
       // call the line is okay
      point_line_gradient(P, A, B, gradient);
    }
  }
  //---------------------------------------------------------------------------

  pair<double, bool> 
  point_line(const std::vector<double>& P, const std::vector<double>& A,
                    const std::vector<double>& B, const std::string type)
  {
    std::vector<double> I;
    return point_line(P, A, B, type, I);
  }
  //---------------------------------------------------------------------------
  
  std::pair<double, bool> 
  point_line(const std::vector<double>& P, const std::vector<double>& A,
             const std::vector<double>& B, const std::string type,
             std::vector<double>& I)
  {
    assert(type == "sd" or type == "d");
    assert(P.size() == A.size() and A.size() == B.size());

    // compute line parameter t of the intersect and the intersect
    const double AB = norm(B - A, 2);

    assert(not compare(AB, 0.0)); // |AB-0.0| < LA_EPS is false 

    const double t = dot(P - A, B - A)/AB/AB;
    bool inside = (t > 1 or t < 0) ? false : true;
    I = A + t*(B - A);
    if(type == "d") 
    {
      return pair<double, double>(norm(P - I, 2), inside);
    }
    else
    {
      // check the sign
      int sign = cross(I-A, P - I)[2] >= 0 ? 1 : -1;
      return pair<double, double>(sign*norm(P - I, 2), inside);
    }
  }
  //---------------------------------------------------------------------------
  
  void point_line_gradient(const std::vector<double>& P,
                           const std::vector<double>& A,
                           const std::vector<double>& B,
                           std::vector<double>& gradient)
  {
    // get the distance and intersect
    std::vector<double> I;
    std::pair<double, bool> distance_onEdge=point_line(P, A, B, "d", I);
    gradient = (P-I)/distance_onEdge.first;
  }
  //---------------------------------------------------------------------------

  double point_polygon(const std::vector<double>& P,
                       const std::vector<double>& vertices,
                       const std::string type)
  {
    assert(type == "sd" or type == "d");
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

    //
    if(type == "d")
    {
      return *min_element(distances.begin(), distances.end());
    }
    else
    {
      const int sign = g_inside(vertices, P) ? -1 : 1;
      return sign*(*min_element(distances.begin(), distances.end()));
    }
  }
  //---------------------------------------------------------------------------
  
  void point_polygon_gradient(const std::vector<double>& P,
                                const std::vector<double>& vertices,
                                std::vector<double>& _gradient)
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
      else
      {
        distances.push_back(1E8); // something big too keep the logic below OK
      }
    }

    // find index if the min element
    double minima = *min_element(distances.begin(), distances.end());
    std::size_t i = std::distance(distances.begin(),
                         std::find(distances.begin(), distances.end(), minima));
    // we have ordered distances as vertex, edge from that vertex, vertex ...
    // if i even it is a vertex
    if((i%2) == 0)
    {
      const double* _Q = &vertices[(i/2)*dim];
      const std::vector<double> Q(_Q, _Q + dim);
      point_point_gradient(P, Q, _gradient);
    }
    else
    {
      const double* _A = &vertices[(i/2)*dim];
      const std::vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[(((i/2)+1)%n_edges)*dim];
      const std::vector<double> B(_B, _B + dim);

      point_edge_gradient(P, A, B, _gradient);
    }
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
  
  std::pair<double, bool> point_plane(const std::vector<double>& P,
                                      const std::vector<double>& A,
                                      const std::vector<double>& B,
                                      const std::vector<double>& C,
                                      const std::string type)
  {
    assert(type == "sd" or type == "d");
    assert(P.size() == 3 and A.size() == 3 and B.size() == 3 and C.size() == 3);

    // set up the matrix problem [[a, b], [c, d]][[t],[s]] = [[e], [f]]
    // to get t,s such that P-I, I = At + Cs + B(1-t-s) is perp to A-B and C-B
    const double a = dot(A - B, A - B);
    const double b = dot(A - B, C - B);
    const double c = b;
    const double d = dot(C - B, C - B);
    const double e = dot(P - B, A - B);
    const double f = dot(P - B, C - B);
    
    const double det = a*d - c*b;
    assert(not compare(det, 0.0));

    const double t = (d*e - b*f)/det;
    const double s = (-c*e + a*f)/det;
    const vector<double> I = A*t + C*s + B*(1 - t - s);
    const bool inside = (0 <= t and t <= 1) and (0 <= s and s <= 1) and
                        (t+s <= 1);
    if(type == "d")
    {
      return pair<double, double>(norm(P - I, 2), inside);
    }
    else
    {
      const int sign = dot(P-I, cross(A-B, C-B)) >= 0 ? 1 : -1;
      return pair<double, double>(sign*norm(P - I, 2), inside);
    }
  }
  //---------------------------------------------------------------------------
  
  double point_3face(const std::vector<double>& P,
                     const std::vector<double>& A,
                     const std::vector<double>& B,
                     const std::vector<double>& C)
  {
    // compute plane params t,s  of the intersect
    const pair<double, bool> d_inside = point_plane(P, A, B, C, "d");
    const bool inside = d_inside.second;
    if(not inside) // intersect outside of plane
    {
      return -1;  // distance should be always positeve, so this is a flag
    }
    else
    {
      return d_inside.first;
    }
  }
  //---------------------------------------------------------------------------
    
  double point_tet(const std::vector<double>& P,
                   const std::vector<double>& vertices,
                   const std::string type)
  {
    assert(type == "sd" or type == "d");
    // minimal distance from the vertices and edges and faces of tet
    const size_t dim = 3;
    assert(P.size() == dim);
    assert(vertices.size()%dim == 0);
    const size_t n_vertices = vertices.size()/dim; // also num faces
    assert(n_vertices == 4);
    const size_t n_edges = 6;

    vector<double> distances; // will hold 2*n_vertices + n_edges values
    distances.reserve(2*n_vertices + n_edges);

    // push back distances from verticess;
    for(size_t i = 0; i < n_vertices; i++)
    {
      const double* _V = &vertices[i*dim];
      const vector<double> V(_V, _V + dim);
      distances.push_back(point_point(P, V));
    }

    // push back distance from faces
    for(size_t i = 0; i < n_vertices; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);

      const double* _B = &vertices[((i+1)%n_vertices)*dim];
      const vector<double> B(_B, _B + dim);

      const double* _C = &vertices[((i+2)%n_vertices)*dim];
      const vector<double> C(_C, _C + dim);
      const double d = point_3face(P, A, B, C);
      if(d != -1)
      {
        distances.push_back(d);
      }
    }

    // push back edge distances
    // first cover vertices by sliding filter: 01, 12, 23, 30
    for(size_t i = 0;  i < n_vertices; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[((i+1)%n_vertices)*dim];
      const vector<double> B(_B, _B + dim);
      const double d = point_edge(P, A, B);
      if(d != -1)
      {
        distances.push_back(d);
      }
    }

    // cover remaining edges: 02 13
    for(size_t i = 0;  i < 2; i++)
    {
      const double* _A = &vertices[i*dim];
      const vector<double> A(_A, _A + dim);
      
      const double* _B = &vertices[(i+2)*dim];
      const vector<double> B(_B, _B + dim);

      const double d = point_edge(P, A, B);
      if(d != -1)
      {
        distances.push_back(d);
      }
    }
    
    if(type == "d")
    {
      return *min_element(distances.begin(), distances.end());
    }
    else
    {
      const int sign = t_inside(vertices, P) ? -1 : 1;
      return sign*(*min_element(distances.begin(), distances.end()));
    }
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
