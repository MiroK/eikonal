#include "RectangleMeshGenerator.h"
#include <dolfin/mesh/Mesh.h>
#include <dolfin/generation/RectangleMesh.h>
#include <cmath>

using namespace dolfin;

namespace eikonal
{
  RectangleMeshGenerator::RectangleMeshGenerator(const std::size_t _i_min,
                                                 const std::size_t _i_max,
                                    const double _LL_x, const double _LL_y,
                                    const double _UR_x, const double _UR_y) :
  i_min(_i_min), i_max(_i_max), i(_i_min), LL_x(_LL_x), LL_y(_LL_y),
  UR_x(_UR_x), UR_y(_UR_y) { }
  //----------------------------------------------------------------------------

  bool RectangleMeshGenerator::end() const { return i == i_max; }
  //----------------------------------------------------------------------------

  boost::shared_ptr<dolfin::Mesh> RectangleMeshGenerator::operator*() const
  {
    std::size_t N = (std::size_t)pow(2, i);
    return boost::shared_ptr<RectangleMesh>
            (new RectangleMesh(LL_x, LL_y, UR_x, UR_y, N, N));
  }
  //----------------------------------------------------------------------------

  void RectangleMeshGenerator::operator++() { i++; }
  //----------------------------------------------------------------------------

  std::string RectangleMeshGenerator::type() const
  { return std::string("f"); } 
  //----------------------------------------------------------------------------
}
