#include "test.h"
#include <dolfin/function/Function.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>

using namespace dolfin;

void LowerBoundarySeeder::initialize(Function& u, MeshFunction<bool>& mesh_f,
                                     std::vector<std::size_t>& fixed_dofs,
                                     std::string type)
{
  // there is no signed distance for this test
  if(type == "sd")
    error("test.cpp", "failed to initialize lower boundary",
          "only distance function can be created");


}

