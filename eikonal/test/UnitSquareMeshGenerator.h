#ifndef _UNIT_SQUARE_MESH_GENERATOR_H_
#define _UNIT_SQUARE_MESH_GENERATOR_H_

/*
  Generate UnitSquareMeshes(2*i, 2*i) with i in [i_min, i_max).
*/

#include "RectangleMeshGenerator.h"

namespace eikonal
{
  class UnitSquareMeshGenerator : public RectangleMeshGenerator
  {
  public:
    // constructor, set the bounds on i
    UnitSquareMeshGenerator(const std::size_t _i_min,
                            const std::size_t _i_max,
                            const bool _perturbed=false);
  };
}


#endif // _UNIT_SQUARE_MESH_GENERATOR_H_
