#ifndef _RECTANGLE_MESH_GENERATOR_H_
#define _RECTANGLE_MESH_GENERATOR_H_

/*
  Generator RectangleMeshes (same syntax as dolfin on box bounds). Number of
  division same in x, y directions = 2^i with i_min <= i < i_max
*/

#include "MeshGenerator.h"

namespace eikonal
{
  class RectangleMeshGenerator : public MeshGenerator
  {
  public:
    // constructor, set the bounds on i
    RectangleMeshGenerator(const std::size_t _i_min,
                           const std::size_t _i_max,
                           const double _LL_x, const double _LL_y,
                           const double _UR_x, const double _UR_y);

    // iteration is over
    virtual bool end() const;

    // get the current mesh
    virtual boost::shared_ptr<dolfin::Mesh> operator*() const;      

    // increment
    virtual void operator++();

    // get the type of mesh "f" for fenics meshes, "g" for meshes from gmsh
    virtual std::string type() const;
  
  private:
    const double LL_x;
    const double LL_y;
    const double UR_x;
    const double UR_y;

    const std::size_t i_min;
    const std::size_t i_max;
    std::size_t i;
  };
}

#endif // _RECTANGLE_MESH_GENERATOR_H_
