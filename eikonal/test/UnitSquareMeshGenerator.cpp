#include "UnitSquareMeshGenerator.h"

using namespace dolfin;

namespace eikonal
{
  UnitSquareMeshGenerator::UnitSquareMeshGenerator
  (const std::size_t _i_min, const std::size_t _i_max, const bool _perturbed) : 
  RectangleMeshGenerator(_i_min, _i_max, 0., 0., 1., 1., _perturbed) { }
  //----------------------------------------------------------------------------
}
