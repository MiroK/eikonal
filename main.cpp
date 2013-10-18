#include "eikonal.h"

using namespace eikonal;

int main()
{
  /* 2 circle 
  double _c1[2] = {-1., 0.}; std::vector<double> c1(_c1, _c1+2);
  double _c2[2] = {sqrt(1.5), 0}; std::vector<double> c2(_c2, _c2+2);
  TwoCircles two_circles("two_circles", c1, 0.5, c2, 0.5);
  Problem problem(two_circles);*/

  /* polygon
  double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
  Polygon polygon("triangle", c, 1, 3);
  Problem problem(polygon);*/

  /* zalesak
  double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
  const double R = 0.5, W = 0.25, L = 0.75;
  Zalesak zalesak(c, R, W, L, 1000);
  Problem problem(zalesak);*/

  // line
  double _A[2] = {-2., -2.}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {2., -2}; std::vector<double> B(_B, _B+2);
  Segment segment("segment", A, B);
  Problem problem(segment);


  dolfin::RectangleMesh mesh(-2, -2, 2, 2, 100, 100);

  int status;
  //perform convergence test on problem using Solver on fenics' UnitSquareMesh
  RectangleMeshGenerator mesh_gen(3, 8, -2, -2, 2, 2, false);

  status = linear_2D_test<Solver>(problem, mesh_gen, true);

  return 0;
}
