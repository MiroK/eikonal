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
  std::vector<la_index> fixed_dofs(_fs, _fs + 4);
  
  Function u(V);
  for(std::size_t i = 0; i < 16; i++)
    u.vector()->setitem(i, 100);

  for(std::size_t i = 0; i < 4; i++)
    u.vector()->setitem(fixed_dofs[i], 0);
 
  plot(u);
  interactive(true);
  
  std::size_t iter_count = solver.solve(u, fixed_dofs);
  std::cout << "Eikonal solver finished in " << iter_count << " iterations.\n";
  
  plot(u);
  interactive(true);*/

  double _A[3] = {0, 0, 0}; std::vector<double> A(_A, _A + 3);
  double _B[3] = {1, 0, 0}; std::vector<double> B(_B, _B + 3);
  Segment segment(A, B);

  std::vector<Point> points;
  segment.seed(points, 4);
  
  return 0;
}
