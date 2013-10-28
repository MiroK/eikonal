#include "GmshMeshGenerator.h"
#include <sstream>
#include <dolfin/mesh/Mesh.h>

using namespace dolfin;

namespace eikonal
{
  GmshMeshGenerator::GmshMeshGenerator(const std::size_t _i_min,
                                       const std::size_t _i_max,
                                       const std::string _root,
                                       const bool _smoothed) :
  i_min(_i_min), i_max(_i_max), i(_i_min), root(_root), smoothed(_smoothed) { }
  //----------------------------------------------------------------------------
  
  bool GmshMeshGenerator::end() const { return i == i_max; }
  //----------------------------------------------------------------------------

  boost::shared_ptr<dolfin::Mesh>
  GmshMeshGenerator::operator*()     
  {
    std::ostringstream help;
    help << "/home/miro3/Documents/Programming/Cpp/Eikonal/test_results/meshes/"
    << root.c_str() << "_" << i << ".msh.xml";
    std::string mesh_name = help.str();
    boost::shared_ptr<Mesh> mesh(new Mesh(mesh_name));

    if(smoothed)
    {
      mesh->smooth();
    }
    return mesh;
  }
  //----------------------------------------------------------------------------

  void GmshMeshGenerator::operator++() { i++; } 
  //----------------------------------------------------------------------------

  std::string GmshMeshGenerator::type() const
  { return smoothed ? std::string("gsmooth") : std::string("g"); }
  //----------------------------------------------------------------------------
}
