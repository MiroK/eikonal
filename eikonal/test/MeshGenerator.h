#ifndef _MESH_GENERATOR_H_
#define _MESH_GENERATOR_H_

/*
  Create meshes on the fly.
*/

#include <boost/shared_ptr.hpp>
#include <string>

namespace dolfin
{
  class Mesh;
}

namespace eikonal
{
  class MeshGenerator
  {
  public:
    // iteration is over
    virtual bool end() const = 0;

    // get the current mesh
    virtual boost::shared_ptr<dolfin::Mesh> operator*() = 0;      

    // increment
    virtual void operator++() = 0;

    // get the type of mesh "f" for fenics meshes, "g" for meshes from gmsh
    virtual std::string type() const = 0;
  };
}

#endif // _MESH_GENERATOR_H_
