#include "eikonal.h"
#include <cstdlib>
#include "test/CG1.h"
#include <dolfin.h>

using namespace eikonal;
using namespace dolfin;

int main(int argc, char* argv[])
{
  /*double S[2] = {0.25, -0.5};
  local_test_S(S);*/

  UnitSquareMesh mesh(2, 2);
  CG1::FunctionSpace V(mesh);
  t_smap_la cell_dof = cell_to_dof(V);
  la_vmap_t dof_cell = dof_to_cell(cell_dof);
  la_smap_la dof_dof = dof_to_dof(cell_dof, dof_cell);

  std::cout << print_map<la_smap_la, std::set<dolfin::la_index> >(dof_dof, "x");

/*
  Dolphin dolphin;
  Problem problem(dolphin);*/
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

