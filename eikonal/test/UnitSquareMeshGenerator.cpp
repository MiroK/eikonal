#include "UnitSquareMeshGenerator.h"
#include <dolfin/mesh/Mesh.h>
#include <dolfin/generation/UnitSquareMesh.h>
#include <cmath>

using namespace dolfin;

namespace eikonal
{
  UnitSquareMeshGenerator::UnitSquareMeshGenerator
  (const std::size_t _i_min, const std::size_t _i_max) : i_min(_i_min),
  i_max(_i_max), i(_i_min) { }
  //----------------------------------------------------------------------------

  bool UnitSquareMeshGenerator::end() const { return i == i_max; }
  //----------------------------------------------------------------------------

  boost::shared_ptr<dolfin::Mesh>
  UnitSquareMeshGenerator::operator*() const     
  {
    std::size_t N = (std::size_t)pow(2, i);
    return boost::shared_ptr<UnitSquareMesh>(new UnitSquareMesh(N, N));
  }
  //----------------------------------------------------------------------------

  void UnitSquareMeshGenerator::operator++() { i++; } 
  //----------------------------------------------------------------------------

  std::string UnitSquareMeshGenerator::type() const { return std::string("f"); }
  //----------------------------------------------------------------------------
}
