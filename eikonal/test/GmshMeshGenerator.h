#ifndef _GMSH_MESH_GENERATOR_H_
#define _GMSH_MESH_GENERATOR_H_

/*
  Generate meshes from .msh files.
*/

#include "MeshGenerator.h"
#include <string>

namespace eikonal
{
  class GmshMeshGenerator : public MeshGenerator
  {
  public:
    // constructor, set the bounds on i and root
    // names are constructed as root_i.msh.xml
    GmshMeshGenerator(const std::size_t _i_min,
                      const std::size_t _i_max,
                      const std::string _root);

    // iteration is over
    virtual bool end() const;

    // get the current mesh
    virtual boost::shared_ptr<dolfin::Mesh> operator*();      

    // increment
    virtual void operator++();

    // get the type of mesh "f" for fenics meshes, "g" for meshes from gmsh
    virtual std::string type() const;
  
  private:
    const std::size_t i_min;
    const std::size_t i_max;
    const std::string root;
    std::size_t i;
  };
}

#endif // _GMSH_MESH_GENERATOR_H_
