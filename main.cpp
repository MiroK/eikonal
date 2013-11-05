#include "eikonal.h"
#include <cstdlib>
#include "test/CG1.h"
#include <dolfin.h>

using namespace eikonal;
using namespace dolfin;

int main(int argc, char* argv[])
{
  double _P[2] = {1., 0.}; std::vector<double> P(_P, _P+2);
  MyPoint point(P);
  Problem problem(point);
 
  UnitSquareMesh mesh(2, 2);
  CG1::FunctionSpace V(mesh);

  dolfin::Function u(V);
  dolfin::Function u_exact(V);
  std::set<dolfin::la_index> fixed_dofs;
      
  problem.init(fixed_dofs, u);
  problem.exact_solution(u_exact);

  FMMSolver solver(V);
  solver.solve(u, fixed_dofs); 
  
  plot(u);
  interactive(true);

  //int status;
  // convergence test on meshes by gmsh 0 .. 6, smoothing
  //GmshMeshGenerator mesh_gen3(6, 7, "rectangle", true);
  //status = hermite_test(problem, mesh_gen3, 2, "corners", 2, false);
  
  //assert(argc == 3);
  // LINEAR TESTS
  //all_linear_tests(atoi(argv[1]), atoi(argv[2]));
 
  // HERMITE TESTS
  //enum HERTMITE_SOLVER {CORNERS, SURFACE, DISTANCE};
  /*std::size_t solver_choice = atoi(argv[1]);
  std::size_t p_norm = atoi(argv[2]);

  if(solver_choice == CORNERS)
  {
    all_hermite_tests("corners", p_norm);  
  }*/
  /*else if(solver_choice == SURFACE)
  {
    all_hermite_tests("surface", p_norm);
  }
  else if(solver_choice == DISTANCE)
  {
    all_hermite_tests("distance", p_norm);
  }
  */
  return 0;
}

