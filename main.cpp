#include "eikonal.h"
#include <cstdlib>

using namespace eikonal;

int main(int argc, char* argv[])
{
  double S[2] = {0.25, -0.5};
  local_test_S(S);



/*
  Dolphin dolphin;
  Problem problem(dolphin);
  int status;
  // convergence test on meshes by gmsh 0 .. 6, smoothing
  //GmshMeshGenerator mesh_gen3(6, 7, "rectangle", true);
  //status = hermite_test(problem, mesh_gen3, 2, "corners", 2, false);
  
  assert(argc == 3);
  // LINEAR TESTS
  //all_linear_tests(atoi(argv[1]), atoi(argv[2]));
 
  // HERMITE TESTS
  enum HERTMITE_SOLVER {CORNERS, SURFACE, DISTANCE};
  std::size_t solver_choice = atoi(argv[1]);
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

