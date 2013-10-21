#include "linear_test.h"
#include "Seeder.h"
#include "Problem.h"
#include "GmshMeshGenerator.h"
#include "RectangleMeshGenerator.h"
#include "gs/Solver.h"
#include "gs/LinMinSolver.h"
#include "gs/LinNewtonSolver.h"

namespace eikonal
{
  int all_linear_tests(int solver, std::size_t precision)
  {
    enum SOLVERS { LIN_GEOMETRIC, LIN_BRENT, LIN_NEWTON};

    if(solver == LIN_GEOMETRIC)
    {
      std::cout << "Solving with linear geometric solver:" << std::endl;
      run_linear_test<Solver>("line");
      run_linear_test<Solver>("point");
      run_linear_test<Solver>("twocircle");
      run_linear_test<Solver>("triangle");
      run_linear_test<Solver>("zalesak");
      run_linear_test<Solver>("dolphin");
    }
    if(solver == LIN_BRENT)
    {
      std::cout << "Solving with linear Brent solver" << std::endl;
      std::cout << precision << std::endl;
      run_linear_test<LinMinSolver>("line", precision);
      run_linear_test<LinMinSolver>("point", precision);
      run_linear_test<LinMinSolver>("twocircle", precision);
      run_linear_test<LinMinSolver>("triangle", precision);
      run_linear_test<LinMinSolver>("zalesak", precision);
      run_linear_test<LinMinSolver>("dolphin", precision);
    }

    if(solver == LIN_NEWTON)
    {
      std::cout << "Solving with linear Newton solver:" << std::endl;
      std::cout << precision << std::endl;
      run_linear_test<LinNewtonSolver>("line", precision);
      run_linear_test<LinNewtonSolver>("point", precision);
      run_linear_test<LinNewtonSolver>("twocircle", precision);
      run_linear_test<LinNewtonSolver>("triangle", precision);
      run_linear_test<LinNewtonSolver>("zalesak", precision);
      run_linear_test<LinNewtonSolver>("dolphin", precision);
    }
    return 1;
  }
  //---------------------------------------------------------------------------
}
