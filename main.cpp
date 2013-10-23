#include "eikonal.h"
#include <cstdlib>

using namespace eikonal;

int main(int argc, char* argv[])
{

  double _P[2] = {1.0, 0.0}; std::vector<double> P(_P, _P+2);
  MyPoint point(P);
  Problem problem(point);
  
  int status;
  std::size_t precision = 2;
  std::string ordering("distance");
  std::size_t p = 2;
  // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
  UnitSquareMeshGenerator mesh_gen0(1, 4, false);
  status = hermite_test(problem, mesh_gen0, precision, ordering, p, false);
  //std::cout << "\nlinear\n";
  //UnitSquareMeshGenerator mesh_gen1(1, 4, false);
  //status = linear_2D_test<Solver>(problem, mesh_gen1, precision, false);
  //assert(argc == 3);
  // LINEAR TESTS
  //all_linear_tests(atoi(argv[1]), atoi(argv[2]));
 
  //assert(argc == 3);
  // HERMITE TESTS
  /*enum HERTMITE_SOLVER {CORNERS, SURFACE, DISTANCE};
  std::size_t solver_choice = atoi(argv[1]);
  std::size_t p_norm = atoi(argv[2]);

  if(solver_choice == CORNERS)
  {
    all_hermite_tests("corners", p_norm);  
  }
  else if(solver_choice == SURFACE)
  {
    all_hermite_tests("surface", p_norm);
  }
  else if(solver_choice == DISTANCE)
  {
    all_hermite_tests("distance", p_norm);
  }*/
  
  return 0;
}

