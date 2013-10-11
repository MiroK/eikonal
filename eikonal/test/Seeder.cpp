#include "Seeder.h"
#include "la/la_loop.h"
#include "cg/cg_distance.h"
#include <dolfin/mesh/Point.h>
#include <cassert>
#include <cmath>

#include <dolfin/log/dolfin_log.h>
#include "la/la_common.h"

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
}
