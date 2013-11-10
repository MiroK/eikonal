#include "linear_test.h"
#include "Seeder.h"
#include "Problem.h"
#include "GmshMeshGenerator.h"
#include "RectangleMeshGenerator.h"
#include "gs/Solver.h"
#include "gs/LinMinSolver.h"
#include "gs/LinNewtonSolver.h"
#include "gs/FMMSolver.h"

namespace eikonal
{
  int all_linear_tests(int solver, std::size_t precision)
  {
    enum SOLVERS {LIN_GEOMETRIC, LIN_BRENT, LIN_NEWTON, FMM_GEOMETRIC};

    if(solver == LIN_GEOMETRIC)
    {
      std::cout << "Solving with linear geometric solver:" << std::endl;
      run_linear_test<Solver>("point");
      run_linear_test<Solver>("twocircle");
      run_linear_test<Solver>("triangle");
      run_linear_test<Solver>("zalesak");
    }
    if(solver == LIN_BRENT)
    {
      std::cout << "Solving with linear Brent solver" << std::endl;
      std::cout << precision << std::endl;
      run_linear_test<LinMinSolver>("point", precision);
      run_linear_test<LinMinSolver>("twocircle", precision);
      run_linear_test<LinMinSolver>("triangle", precision);
      run_linear_test<LinMinSolver>("zalesak", precision);
    }

    if(solver == LIN_NEWTON)
    {
      std::cout << "Solving with linear Newton solver:" << std::endl;
      std::cout << precision << std::endl;
      run_linear_test<LinNewtonSolver>("point", precision);
      run_linear_test<LinNewtonSolver>("twocircle", precision);
      run_linear_test<LinNewtonSolver>("triangle", precision);
      run_linear_test<LinNewtonSolver>("zalesak", precision);
    }

    if(solver == FMM_GEOMETRIC)
    {
      std::cout << "Solving with FMM geometric:" << std::endl;
      std::cout << precision << std::endl;
      run_linear_test<FMMSolver>("point", precision);
      run_linear_test<FMMSolver>("twocircle", precision);
      run_linear_test<FMMSolver>("triangle", precision);
      run_linear_test<FMMSolver>("zalesak", precision);
    }
    return 1;
  }
  //---------------------------------------------------------------------------
}
