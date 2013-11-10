#include "hermite_test.h"
#include "gs/HermiteSolver.h"

namespace eikonal
{
  int all_hermite_tests(std::size_t solver_type,
                        std::string ordering, std::size_t p)
  {
    enum SOLVERS {FMM, FSM};
 
    if(solver_type == FSM)
    {
      std::cout << "Solving with fsm hermite solver:" << std::endl;
      std::cout << "ordering = " << ordering << " L^p = " << p << std::endl;

      std::size_t precision = 2;
      run_hermite_test<HermiteSolver>("point", precision, ordering, p); 
      /*run_hermite_test<HermiteSolver>("twocircle", precision, ordering, p); 
      run_hermite_test<HermiteSolver>("triangle", precision, ordering, p); 
      run_hermite_test<HermiteSolver>("zalesak", precision, ordering, p);*/ 
    }    
    return 0;
  }
  //---------------------------------------------------------------------------

}
