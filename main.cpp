#include "eikonal.h"
#include <dolfin.h>

using namespace eikonal;

int main()
{
  // define the seeder should be [0, 0] and [1, 0] but eps must be added
  // in y-coord to help numerics
  const double eps = .001;
  double _A[2] = {0, eps}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1, eps}; std::vector<double> B(_B, _B + 2);
  print(A);
  print(B);
  Segment segment("segment_[0,0]_[1,0]", A, B);

  // define the problem
  //Problem problem(segment);

  /* test 2 circle problem
  double _c1[2] = {-1., -1.}; std::vector<double> c1(_c1, _c1+2);
  double _c2[2] = {1., 1.}; std::vector<double> c2(_c2, _c2+2);
  TwoCircles two_circles("two_circles", c1, 0.5, c2, 0.5);

  
  Problem problem(two_circles);

  dolfin::RectangleMesh mesh(-2, -2, 2, 2, 100, 100);
  CG1::FunctionSpace V(mesh);
  dolfin::Function u(V);
  problem.exact_solution(u);
  dolfin::plot(u);
  dolfin::interactive(true); */

  
  UnitSquareMeshGenerator x(2, 3, true); // 2**4 to 2**4, box [0,0]x[1,1]
  dolfin::plot(*(*x));
  dolfin::interactive(true);
  print((*x)->coordinates());
  

  // int status;
  // perform convergence test on problem using Solver on fenics' UnitSquareMesh
  // UnitSquareMeshGenerator sq_mesh_gen(3, 5);
  // status = linear_2D_test<Solver>(problem, sq_mesh_gen, false);

  // perform convergence test on problem using Solver on gmsh meshes
  // GmshMeshGenerator g_mesh_gen(0, 5, "sqr");
  // status = linear_2D_test<Solver>(problem, g_mesh_gen, true);

  return 0;
}
