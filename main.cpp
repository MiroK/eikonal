#include "eikonal.h"

using namespace eikonal;

int main()
{
  /* test 2 circle problem
  double _c1[2] = {-1., 0.}; std::vector<double> c1(_c1, _c1+2);
  double _c2[2] = {sqrt(1.5), 0}; std::vector<double> c2(_c2, _c2+2);
  TwoCircles two_circles("two_circles", c1, 0.5, c2, 0.5);
  Problem problem(two_circles);*/

  double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
  Polygon polygon("triangle", c, 1, 3);
  Problem problem(polygon);



  dolfin::RectangleMesh mesh(-2, -2, 2, 2, 100, 100);

  int status;
  //perform convergence test on problem using Solver on fenics' UnitSquareMesh
  RectangleMeshGenerator mesh_gen(3, 8, -2, -2, 2, 2, false);

  status = linear_2D_test<Solver>(problem, mesh_gen, false);

  return 0;
}
