#include "GmshMeshGenerator.h"
#include <sstream>
#include <dolfin/mesh/Mesh.h>

using namespace dolfin;

namespace eikonal
{
  GmshMeshGenerator::GmshMeshGenerator(const std::size_t _i_min,
                                       const std::size_t _i_max,
                                       const std::string _root) :
  i_min(_i_min), i_max(_i_max), i(_i_min), root(_root) { }
  //----------------------------------------------------------------------------
  
  bool GmshMeshGenerator::end() const { return i == i_max; }
  //----------------------------------------------------------------------------

  boost::shared_ptr<dolfin::Mesh>
  GmshMeshGenerator::operator*()     
  {
    std::ostringstream help;
    help << root.c_str() << "_" << i << ".msh.xml";
    std::string mesh_name = help.str();
    return boost::shared_ptr<Mesh>(new Mesh(mesh_name));
  }
  //----------------------------------------------------------------------------

  void GmshMeshGenerator::operator++() { i++; } 
  //----------------------------------------------------------------------------

  std::string GmshMeshGenerator::type() const { return std::string("g"); }
  //----------------------------------------------------------------------------
}
