#include "CG1.h"
#include "CG2.h"
#include <dolfin/generation/UnitSquareMesh.h>
#include <iostream>

#include "connectivity.h"
#include "test.h"
#include "linear_algebra.h"
#include "polyhedron.h"

using namespace dolfin;

int main()
{
  /*UnitSquareMesh mesh(100, 100);
  CG2::FunctionSpace V(mesh);
  
  Function u(V);
  MeshFunction<bool> mesh_f(mesh, 2);
  std::vector<std::size_t> fixed_dofs;

  LowerBoundarySeeder lb_seeder;
  lb_seeder.initialize(u, mesh_f, fixed_dofs, "d");

   build all mappings
  t_smap_la _cell_to_dof = cell_to_dof(V);
  la_vmap_t _dof_to_cell = dof_to_cell(_cell_to_dof);
  t_smap_t _cell_to_facet = cell_to_facet(V);
  t_vmap_t _facet_to_cell = facet_to_cell(_cell_to_facet);
  t_smap_t _facet_to_dof = facet_to_dof(V);
  t_vmap_t _dof_to_facet = dof_to_facet(_facet_to_dof);
  t_vmap_d _dof_to_coordinate = dof_to_coordinate(V);
  */

  double _u[3] = {1., 0., 0.};
  double _v[3] = {0., 1, 0.};
  std::vector<double> u(_u, _u + 3);
  std::vector<double> v(_v, _v + 3);
  
  // test dot
  double a = dot(u, v);
  info("%g", a);

  // test cross
  std::vector<double> w = cross(u, v);
  info("%g %g %g", w[0], w[1], w[2]);

  std::vector<double> uvw;
  uvw.insert(uvw.begin(), u.begin(), u.end());
  uvw.insert(uvw.begin() + 3, v.begin(), v.end());
  uvw.insert(uvw.begin() + 6, w.begin(), w.end());
  info("size %d", uvw.size());
  std::vector<double> uvw_c = barycenter(uvw, 3);
  info("%g %g %g", uvw_c[0], uvw_c[1], uvw_c[2]);
  return 0;
}
