#include "Seeder.h"
#include "la/la_loop.h"
#include "cg/cg_distance.h"
#include <dolfin/mesh/Point.h>
#include <cassert>

#include <dolfin/log/dolfin_log.h>
#include "la/la_common.h"

using namespace dolfin;

namespace eikonal
{
  Segment::Segment(const std::vector<double>& _A, const std::vector<double>& _B)
                          : A(_A), B(_B), dim(_A.size())
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
}
