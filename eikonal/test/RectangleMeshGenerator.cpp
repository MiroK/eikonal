#include "RectangleMeshGenerator.h"
#include <dolfin/mesh/Mesh.h>
#include <dolfin/generation/RectangleMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/DomainBoundary.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>

using namespace dolfin;

namespace eikonal
{
RectangleMeshGenerator::RectangleMeshGenerator(const std::size_t _i_min,
                                               const std::size_t _i_max,
                                  const double _LL_x, const double _LL_y,
                                  const double _UR_x, const double _UR_y,
                                  const bool _perturbed) :
i_min(_i_min), i_max(_i_max), i(_i_min), LL_x(_LL_x), LL_y(_LL_y),
UR_x(_UR_x), UR_y(_UR_y), perturbed(_perturbed), deref_count(0)
{
  assert(i_min < i_max);
}
//----------------------------------------------------------------------------

bool RectangleMeshGenerator::end() const { return i == i_max; }
//----------------------------------------------------------------------------

boost::shared_ptr<dolfin::Mesh> RectangleMeshGenerator::operator*()
{
  assert(i < i_max);
  deref_count++;

  if((deref_count-1) > (i-i_min))
  { // mesh is already created, just return it
    return mesh;
  }

  // craete new mesh 
  std::size_t N = (std::size_t)pow(2, i);
  mesh = boost::shared_ptr<RectangleMesh>
           (new RectangleMesh(LL_x, LL_y, UR_x, UR_y, N, N, "crossed"));
  if(not perturbed)
  { // no perturabations needed
    return mesh;
  }
  else
  {
    // every internal vertex should be shifted by a random vector
    // mark facets on the boundary true
    DomainBoundary domain_boundary;
    MeshFunction<bool> boundary_facets(*mesh, mesh->topology().dim()-1);
    boundary_facets.set_all(false);
    domain_boundary.mark(boundary_facets, true);

    // mark vertices on the boundary true
    MeshFunction<bool> boundary_vertices(*mesh, mesh->topology().dim()-2);
    boundary_vertices.set_all(false);
    
    for(FacetIterator facet(*mesh); !facet.end(); ++facet)
    {
      if(boundary_facets[*facet])
      {
        for(VertexIterator vertex(*facet); !vertex.end(); ++vertex)
        {
          boundary_vertices[*vertex] = true;
        }
      }
    }
  
    // shift the internal vertices
    srand(time(NULL));
    const double max_shift = 0.125*mesh->hmin();
    std::vector<double>& mesh_coordinates = mesh->coordinates();
    for(VertexIterator vertex(*mesh); !vertex.end(); ++vertex)
    {
      if(not boundary_vertices[*vertex])
      {
        int factor_0 = rand() % 1001;
        int factor_1 = rand() % 1001;
        int sign_0 = (factor_0 % 2) ? 1 : -1;  
        int sign_1 = (factor_1 % 2) ? 1 : -1;  

        const std::size_t i = vertex->index();
        mesh_coordinates[2*i + 0] += sign_0*factor_0*max_shift/1000;
        mesh_coordinates[2*i + 1] += sign_1*factor_1*max_shift/1000;
      }
    }
    return mesh;
  }
}
//----------------------------------------------------------------------------

void RectangleMeshGenerator::operator++() { i++; }
//----------------------------------------------------------------------------

std::string RectangleMeshGenerator::type() const
  { return perturbed ? std::string("fperturbed") : std::string("f"); } 
  //----------------------------------------------------------------------------
}
