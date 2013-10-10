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
  Problem problem(segment);

  int status;
  // perform convergence test on problem using Solver on fenics' UnitSquareMesh
  // UnitSquareMeshGenerator sq_mesh_gen(3, 5);
  // status = linear_2D_test<Solver>(problem, sq_mesh_gen, false);

  // perform convergence test on problem using Solver on gmsh meshes
  GmshMeshGenerator g_mesh_gen(0, 5, "sqr");
  status = linear_2D_test<Solver>(problem, g_mesh_gen, true);

  return 0;
}
