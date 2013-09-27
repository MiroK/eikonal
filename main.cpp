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
  /*MeshFunction<bool> mesh_f(mesh, 2);
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
  
  double _tet[12] = {0.25, 0.25, 0.25,
                     0.5, 0.75, 0.5,
                     0.75, 0.5, 0.5,
                     0.5, 0.5, 0.75};
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
  
  double _A[2] = {0., 0.}; std::vector<double> A(_A, _A+2);
  double _B[2] = {1., 0.}; std::vector<double> B(_B, _B+2);
  double _C[2] = {0., 1.}; std::vector<double> C(_C, _C+2);
  double _D[2] = {1., 1.}; std::vector<double> D(_D, _D+2);
  double _E[2] = {2., -1.}; std::vector<double> E(_E, _E+2);

  print(center);
  std::cout << "distance |AD| (sqrt(2)/2): " << point_point(A, D) << std::endl;
  std::cout << "distance (0,0, 1) " << point_circle(D, center, 1, "d") << std::endl;
  std::cout << "distance " << point_circle(center, center, 1, "d") << std::endl;
  std::cout << "distance " << point_circle(center, center, 1, "sd") << std::endl;
  std::cout << "distance " << point_line(D, B, C, "sd").first << std::endl;
  std::cout << "distance " << point_line(center, B, C, "sd").first << std::endl;
  std::cout << "distance " << point_edge(E, B, C) << std::endl;

  std::vector<double> nul(3);
  nul[0] = 1; nul[1] = 1; nul[2] = 1;
  std::cout << "tet distance" << point_tet(t_barycenter(tet), tet, "sd") << std::endl;
 */ 
  /* test polygon distance
  // create a polygon at 0.5, 0.5, 0.25
  std::size_t dim = 2;
  std::size_t n_vertices = 8;
  std::vector<double> p_center(dim);
  p_center[0]= 0.5; p_center[1] = 0.5;
  std::vector<double> p_vertices = g_vertices(n_vertices, p_center, 0.25);

  UnitSquareMesh mesh(100, 100);
  CG2::FunctionSpace V(mesh);
  Function u(V);
  boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();
  std::vector<double> coords = dofmap->tabulate_all_coordinates(mesh);
  
  size_t n_dofs = V.dim();
  std::vector<double> values(n_dofs); 
  for(std::size_t i = 0; i < n_dofs; i++)
  {
    const std::vector<double> point(&coords[i*dim], &coords[i*dim] + dim);
    values[i] = point_polygon(point, p_vertices, "sd");
  }
  u.vector()->set_local(values);
  plot(u);
  interactive(true);

  // test plane distance
  double _K[3] = {0, 0, 0}; std::vector<double> K(_K, _K + 3);
  double _L[3] = {1, 0, 0}; std::vector<double> L(_L, _L + 3);
  double _M[3] = {0, 1, 0}; std::vector<double> M(_M, _M + 3);
  double _X[3] = {1, 1, 10}; std::vector<double> X(_X, _X + 3);
  double _X1[3] = {-1, -1, 10}; std::vector<double> X1(_X1, _X1 + 3);
  std::cout << point_3face(X, K, L, M) << std::endl;
  std::cout << point_3face(X1, K, L, M) << std::endl;
 
  std::size_t dim = 3;
  UnitCubeMesh mesh(50, 50, 50);
  TET_CG1::FunctionSpace V(mesh);
  Function u(V);
  boost::shared_ptr<const GenericDofMap> dofmap = V.dofmap();
  std::vector<double> coords = dofmap->tabulate_all_coordinates(mesh);
  
  size_t n_dofs = V.dim();
  std::vector<double> values(n_dofs); 
  for(std::size_t i = 0; i < n_dofs; i++)
  {
    const std::vector<double> point(&coords[i*dim], &coords[i*dim] + dim);
    values[i] = point_tet(point, tet, "sd");
  }
  u.vector()->set_local(values);
  
  File file("tet.pvd");
  file << u;
*/

  // local solver comparisons
  // move vertex around to fail on alpha, beta. what happens with minim
  double _C[2] = {0.5, 1}; std::vector<double> C(_C, _C + 2);
  double u_C = 100;

  double _A[2] = {0, 0}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1, 0.}; std::vector<double> B(_B, _B + 2);
  std::vector<double> AB(4);
  AB[0] = A[0]; AB[1] = A[1]; AB[2] = B[0]; AB[3] = B[1];
  
  double _u_AB[2] = {0, 0}; std::vector<double> u_AB(_u_AB, _u_AB + 2);

  std::cout << "geometric " << ls_geometric_2d(C, u_C, AB, u_AB) << std::endl;
  std::cout << "minimize " << ls_minimize_2d(C, u_C, AB, u_AB) << std::endl;
  
  //  move source around to fail on theta, again what happend with minim. method
  double _s[2] = {0.5, -0.5};
  std::vector<double> source(_s, _s + 2);
  u_AB[0] = point_point(A, source);
  u_AB[1] = point_point(B, source);
  
  std::cout << "geometric " << ls_geometric_2d(C, u_C, AB, u_AB) << std::endl;
  std::cout << "minimize " << ls_minimize_2d(C, u_C, AB, u_AB) << std::endl;
  
  return 0;
}
