#include "eikonal.h"
#include <cstdlib>
#include <dolfin.h>
#include "eikonal/test/CG1_VECTOR.h"

using namespace eikonal;

int main(int argc, char* argv[])
{
  dolfin::RectangleMesh mesh(-2, -2, 2, 2, 100, 100);
  CG1::FunctionSpace V(mesh);
 
  dolfin::Function u(V);
  dolfin::Function du_dx(V);
  dolfin::Function du_dy(V);
  
  dolfin::Function exact_u(V);
  dolfin::Function exact_du_dx(V);
  dolfin::Function exact_du_dy(V);
 
  std::set<dolfin::la_index> fixed_dofs;

  /*double _P[2] = {0.5, 0.5}; std::vector<double> P(_P, _P+2);
  MyPoint point(P);
  Problem problem(point);*/
 
  /*
  double _A[2] = {0., 0.}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1., 0.}; std::vector<double> B(_B, _B+2);
  Segment segment("segment", A, B);
  Problem problem(segment);*/
  
  /*
  double _c1[2] = {0.25, 0.25}; std::vector<double> c1(_c1, _c1+2);
  double _c2[2] = {0.75, 0.75}; std::vector<double> c2(_c2, _c2+2);
  TwoCircles two_circles("twocircle", c1, 0.125, c2, 0.125);
  Problem problem(two_circles);*/
  
  /*double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
  Polygon polygon("triangle", c, 1, 3);
  Problem problem(polygon);*/

  /*
  double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
  const double R = 1, W = 0.25, L = 1.5;
  Zalesak zalesak(c, R, W, L, 1000);
  Problem problem(zalesak);*/
  
  Dolphin dolphin;
  Problem problem(dolphin);
  
  problem.init(fixed_dofs, u, du_dx, du_dy);
  problem.exact_solution(exact_u, exact_du_dx, exact_du_dy);

  dolfin::plot(u);
  dolfin::plot(exact_u);
  dolfin::interactive(true);
  /* 
  HermiteSolver solver(V);
  std::size_t precision = 2;
  solver.solve(u, du_dx, du_dy, fixed_dofs, precision);

  dolfin::plot(u);
  dolfin::interactive(true);
  */
  CG1_VECTOR::FunctionSpace VV(mesh);
  dolfin::Function du(VV);
  std::vector<double> du_values(VV.dim());
  
  boost::shared_ptr<dolfin::GenericVector> dx_vector = exact_du_dx.vector();
  boost::shared_ptr<dolfin::GenericVector> dy_vector = exact_du_dy.vector();
  
  
  for(std::size_t i = 0; i < du_values.size()/2; i++)
  {
    du_values[2*i] = (*dx_vector)[i];
    du_values[2*i+1] = (*dy_vector)[i];
  }
  du.vector()->set_local(du_values);
  
  dolfin::plot(du);
  dolfin::interactive(true);

  /*
  CG1_FORMS::Form_norm1 l1(mesh, u, exact_u);
  CG1_FORMS::Form_norm2 l2(mesh, u, exact_u);
  double l1_norm = dolfin::assemble(l1);
  double l2_norm = dolfin::assemble(l2);
  std::cout << mesh.hmin() << " " << l1_norm << " " << l2_norm << std::endl;
  */
  return 0;
}

