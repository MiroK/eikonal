#include <dolfin/generation/UnitSquareMesh.h>
#include <dolfin/generation/RectangleMesh.h>
#include <dolfin/plot/plot.h>
#include "eikonal.h"

using namespace dolfin;
using namespace eikonal;

int main()
{
  /*UnitSquareMesh mesh(3, 3);
  CG1::FunctionSpace V(mesh);

  Solver solver(V);

  la_index _fs[4] = {0, 5, 2, 9};
  std::set<la_index> fixed_dofs(_fs, _fs + 4);
  
  Function u(V);
  for(std::size_t i = 0; i < 16; i++)
    u.vector()->setitem(i, 100);

  std::set<la_index>::const_iterator fixed_dof = fixed_dofs.begin();
  for( ; fixed_dof != fixed_dofs.end(); fixed_dof++)
    u.vector()->setitem(*fixed_dof, 0);
 
  plot(u);
  interactive(true);
  
  std::size_t iter_count = solver.solve(u, fixed_dofs);
  std::cout << "Eikonal solver finished after " << iter_count << " iterations.\n";
  
  plot(u);
  interactive(true);*/             // TEST SWEEP

  /*
  double _A[2] = {0, 0}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1, 0}; std::vector<double> B(_B, _B + 2);
  Segment segment(A, B);

  std::vector<Point> points;
  segment.seed(points, 4);
  for(std::size_t i = 0; i < 4; i++)
    info(points[i].str());

  for(std::size_t i = 0; i < 3; i++)
  {
    for(std::size_t j = 0; j < 3; j++)
    {
      double _point[2];
      _point[0] = 0.5*i;
      _point[1] = 0.5*j;
      std::vector<double> point(_point, _point + 2);
      double distance = segment.distance(point);
      info("distance to [%g,%g] is %g", _point[0], _point[1], distance);
    }
  }*/                               // TEST SEEDER
 
  double _A[2] = {0, 0}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1, 0}; std::vector<double> B(_B, _B + 2);
  Segment segment(A, B);
  Problem problem(segment);
  
  UnitSquareMesh mesh(3, 3);
  CG1::FunctionSpace V(mesh);

  Function u(V);
  std::set<la_index> fixed_dofs;

  problem.init(fixed_dofs, u);

  std::cout << "Fixed dofs(" << fixed_dofs.size() << ") : "; 
  std::set<la_index>::const_iterator fixed_dof = fixed_dofs.begin();
  for(; fixed_dof != fixed_dofs.end(); fixed_dof++) 
    std::cout << *fixed_dof << " ";
  std::cout << std::endl;

  std::cout << "U values: ";
  for(std::size_t i = 0; i < V.dim(); i++)
    std::cout << (*u.vector())[i] << " " << std::endl;
  std::cout << std::endl;

  Function u_exact(V);
  problem.exact_solution(u_exact);

  plot(u_exact);
  interactive(true);
                                       // TEST PROBLEM
  return 0;
}
