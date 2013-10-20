#include "Seeder.h"
#include "la/la_loop.h"
#include "cg/cg_distance.h"
#include "cg/cg_polygon.h"
#include <dolfin/mesh/Point.h>
#include <cassert>
#include <cmath>
#include <fstream>

using namespace dolfin;

namespace eikonal
{
  Segment::Segment(const std::string& name, const std::vector<double>& _A,
                   const std::vector<double>& _B)
                          : Seeder(name), A(_A), B(_B), dim(_A.size())
  {
    assert((_A.size() == _B.size()) and ((dim == 2) or (dim == 3)));
  }
  //---------------------------------------------------------------------------

  void Segment::seed(std::vector<dolfin::Point>& points,
                     const std::size_t num_points) const
  {
    points.resize(num_points);
    const double dt = 1./(num_points-1);
    for(std::size_t i = 0; i < num_points; i++)
    {
      double t = i*dt;
      const std::vector<double> P = A + t*(B-A);
      points[i] = Point(dim, &P.front());
    }
  }
  //---------------------------------------------------------------------------

  double Segment::distance(const std::vector<double>& point) const
  {
    return point_edge(point, A, B);
  }
  //---------------------------------------------------------------------------

  TwoCircles::TwoCircles(const std::string& name,
                         const std::vector<double>& _c1, const double _r1,
                         const std::vector<double>& _c2, const double _r2) :
  Seeder(name), c1(_c1), r1(_r1), c2(_c2), r2(_r2)
  {
    assert((_c1.size() == _c2.size()) and (_c1.size() == 2));
  }
  //---------------------------------------------------------------------------

  void TwoCircles::seed(std::vector<dolfin::Point>& points,
                        const std::size_t num_points) const
  {
    points.reserve(2*num_points);
    const double d_theta = 2.*M_PI/(num_points-1);
    
    // seed the first circle
    for(std::size_t i = 0; i < num_points; i++)
    {
      double P[2];
      P[0] = c1[0] + r1*cos(i*d_theta);
      P[1] = c1[1] + r1*sin(i*d_theta);
      points.push_back(Point(2, &P[0]));
    }

    // seed the second circle
    for(std::size_t i = 0; i < num_points; i++)
    {
      double P[2];
      P[0] = c2[0] + r2*cos(i*d_theta);
      P[1] = c2[1] + r2*sin(i*d_theta);
      points.push_back(Point(2, &P[0]));
    }
  }
  //---------------------------------------------------------------------------
      
  double TwoCircles::distance(const std::vector<double>& point) const
  {
    return std::min(point_circle(point, c1, r1, "d"),
                    point_circle(point, c2, r2, "d"));
  }
  //--------------------------------------------------------------------------

  Polygon::Polygon(const std::string& name, const std::vector<double>& c,
                   const double r, const std::size_t n) : Seeder(name)
  {
    assert(c.size() == 2);
    vertices = g_vertices(n, c, r);
  }
  //---------------------------------------------------------------------------

  void Polygon::seed(std::vector<dolfin::Point>& points,
            const std::size_t num_points) const
  {
    points.clear();

    std::vector<double> points_x = boundary_points(vertices, num_points);
    std::vector<double>::const_iterator point_x = points_x.begin();
    for( ; point_x != points_x.end(); point_x += 2)
    {
      points.push_back(Point(2, &(*point_x)));
    }
  }
  //---------------------------------------------------------------------------

  double Polygon::distance(const std::vector<double>& point) const
  {
    return point_polygon(point, vertices, "d");
  }
  //---------------------------------------------------------------------------
  
  Zalesak::Zalesak(const std::vector<double>& c,
                   const double R, const double W, const double L, 
                   const std::size_t num_vertices) : Seeder("zalesak")
  {
    const double cx = c[0]; 
    const double cy = c[1];

    const double dl = 1./(num_vertices-1);
    const double theta = asin(W/2/R); // angle between axis of notch and notch
                                      // circle intersect
    const double l = R*cos(theta); // dist center to top of notch
    const double len = W + 2*L + 2*R*(M_PI - theta); // length of the disk
    const double x_0 = cx;
    const double y_0 = cy - l + L;

    vertices.resize(2*num_vertices);
    const double d_len = len/num_vertices;
    for(std::size_t i = 0; i < num_vertices; i++)
    {
      const double d = i*d_len;
      
      if(d <= W/2) 
      {
        vertices[2*i + 0] = x_0 + d;
        vertices[2*i + 1] = y_0;
      }
      else if(d <= W/2 + L)
      {
        vertices[2*i + 0] = x_0 + W/2;
        vertices[2*i + 1] = y_0 - (d - W/2);
      }
      else if(d <= len - (W/2 + L))
      {
        const double phi = theta + (d - (W/2 + L))/R;
        vertices[2*i + 0] = cx + R*sin(phi);
        vertices[2*i + 1] = cy - R*cos(phi);
      }
      else if(d <= len - W/2)
      {
        vertices[2*i + 0] = x_0 - W/2;
        vertices[2*i + 1] = y_0 - (len - d - W/2);
      }
      else
      {
        vertices[2*i + 0] = x_0 - (len - d);
        vertices[2*i + 1] = y_0;
      }
    }
  }
  //---------------------------------------------------------------------------

  void Zalesak::seed(std::vector<dolfin::Point>& points,
                     const std::size_t num_points) const
  {
    points.clear();
    
    std::vector<double> points_x = boundary_points(vertices, 3);
    std::vector<double>::const_iterator point_x = points_x.begin();
    for( ; point_x != points_x.end(); point_x += 2)
    {
      points.push_back(Point(2, &(*point_x)));
    }
  }
  //---------------------------------------------------------------------------

  double Zalesak::distance(const std::vector<double>& point) const
  {
    return point_polygon(point, vertices, "d");
  }
  //---------------------------------------------------------------------------

  Dolphin::Dolphin() : Seeder("dolfin")
  {
    std::ifstream file;
    
    file.open("/home/miro3/Documents/Programming/Cpp/Eikonal/test_results/"
              "meshes/dolfin_file.txt", std::ios::in);
    assert(file);

    for(std::size_t i = 0; i < 126; i++)
    {
      for(std::size_t j = 0; j < 2; j++)
      {
        double x;
        file >> x;
        vertices.push_back(x);
      }
    }
  }
  //---------------------------------------------------------------------------

  void Dolphin::seed(std::vector<dolfin::Point>& points,
                    const std::size_t num_points) const
  {
    points.clear();
    
    std::vector<double> points_x = boundary_points(vertices, 3);
    std::vector<double>::const_iterator point_x = points_x.begin();
    for( ; point_x != points_x.end(); point_x += 2)
    {
      points.push_back(Point(2, &(*point_x)));
    }
  }
  //---------------------------------------------------------------------------

  double Dolphin::distance(const std::vector<double>& point) const
  {
    return point_polygon(point, vertices, "d");
  }
  //---------------------------------------------------------------------------

  MyPoint::MyPoint(const std::vector<double>& _vertex) : 
  Seeder("point"), vertex(_vertex) { }
  //---------------------------------------------------------------------------

  void MyPoint::seed(std::vector<dolfin::Point>& points,
                      const std::size_t num_points) const
  {
    points.clear();
    points.push_back(Point(2, &(*vertex.begin())));
  }
  //---------------------------------------------------------------------------

  double MyPoint::distance(const std::vector<double>& point) const
  {
    return point_point(point, vertex);
  }
  //---------------------------------------------------------------------------

  void MyPoint::gradient(const std::vector<double>& point,
                         std::vector<double>& _gradient) const
  {
    point_point_gradient(point, vertex, _gradient); 
  }
  //---------------------------------------------------------------------------
}
