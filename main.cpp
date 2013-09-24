/*#include "CG1.h"
#include "CG2.h"
#include <dolfin/generation/UnitSquareMesh.h>
#include <iostream>

#include "connectivity.h"
#include "test.h"
#include "linear_algebra.h"
#include "polyhedron.h"
#include "polygon.h"*/

#include "eikonal.h"

#include <iostream>
//using namespace dolfin;

using namespace eikonal;

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

  /* test baryc
  std::vector<double> uvw;
  uvw.insert(uvw.begin(), u.begin(), u.end());
  uvw.insert(uvw.begin() + 3, v.begin(), v.end());
  uvw.insert(uvw.begin() + 6, w.begin(), w.end());
  info("size %d", uvw.size());
  std::vector<double> uvw_c = barycenter(uvw, 3);
  info("%g %g %g", uvw_c[0], uvw_c[1], uvw_c[2]);
  */ 
  
  std::vector<double> center(2);
  std::vector<double> square = g_vertices(3, center, 1);
  print(square); 

  std::vector<double> square_center = g_barycenter(square);
  print(square_center);

  std::cout << "boundary: ";
  print(boundary_points(square, 3));

  std::cout << "area: " << area(square) << std::endl;
  std::cout << "inside center " << g_inside(square, square_center) << std::endl;
  g_inside(square, square_center);
  
  std::vector<double> two_zero;
  two_zero.push_back(10);
  two_zero.push_back(11);
  std::cout << "inside [1, 0] " << g_inside(square, two_zero) << std::endl;
  
  double _tet[12] = {0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2};
  std::vector<double> tet(_tet, _tet + 12);
  std::cout << "tet volume: " << volume(tet) << std::endl;
 
  std::vector<double> t_c = t_barycenter(tet);
  print(t_c);

  std::cout << "center inside: " << t_inside(tet, t_c) << std::endl;

  std::vector<double> oo4;
  oo4.push_back(0); oo4.push_back(0); oo4.push_back(4);
  std::cout << "[0, ,0 4] inside: " << t_inside(tet, oo4) << std::endl;

  std::cout << "edges: ";
  print(edge_points(tet, 3));

  return 0;
}
