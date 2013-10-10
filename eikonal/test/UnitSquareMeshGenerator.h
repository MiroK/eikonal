#ifndef _UNIT_SQUARE_MESH_GENERATOR_H_
#define _UNIT_SQUARE_MESH_GENERATOR_H_

/*
  Generate UnitSquareMeshes(2*i, 2*i) with i in [i_min, i_max).
*/

#include "MeshGenerator.h"

namespace eikonal
{
  class UnitSquareMeshGenerator : public MeshGenerator
  {
  public:
    // constructor, set the bounds on i
    UnitSquareMeshGenerator(const std::size_t _i_min, const std::size_t _i_max);

    // iteration is over
    virtual bool end() const;

    // get the current mesh
    virtual boost::shared_ptr<dolfin::Mesh> operator*() const;      

    // increment
    virtual void operator++();

    // get the type of mesh "f" for fenics meshes, "g" for meshes from gmsh
    virtual std::string type() const;
  
  private:
    const std::size_t i_min;
    const std::size_t i_max;
    std::size_t i;
  };
}


#endif // _UNIT_SQUARE_MESH_GENERATOR_H_
