#ifndef _TEST_H_
#define _TEST_H_

#include "Problem.h"
#include "gs/Solver.h"
#include "CG1.h"
#include "CG1_FORMS.h"
#include <dolfin/generation/UnitSquareMesh.h>
#include <dolfin/fem/assemble.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

/*
  Perform convergence test of the Eikonal solver on different problems.
*/

namespace eikonal
{

  // use CG1 in 2D to test Eikonal solvers
  template<typename T> int linear_2D_test(const Problem& problem);
}

// implementation

namespace eikonal
{
  template<typename T> int linear_2D_test(const Problem& problem)
  {
    // 8, 16, 32, 64, 128, 256 <--> 3, 4, 5, 6, 7, 8
    std::ofstream file;
    std::string file_name;
    file.open(file_name.c_str(), std::ios::out);
    file.close();
    for(std::size_t i = 3; i < 5; i++)
    {
      file.open(file_name.c_str(), std::ios::app);

      std::size_t N = (std::size_t)pow(2, i);
      dolfin::UnitSquareMesh mesh(N, N);
      CG1::FunctionSpace V(mesh);
      
      dolfin::Function u(V);
      dolfin::Function u_exact(V);
      std::set<dolfin::la_index> fixed_dofs;
      
      // set boundary condition, u, and get the exact solution
      problem.init(fixed_dofs, u);
      problem.exact_solution(u_exact);

      // create solver and compute the solution
      T solver(V);
      std::size_t num_iters = solver.solve(u, fixed_dofs); 
  
      // get the error of in L1, L2 norms
      CG1_FORMS::Form_norm1 l1(mesh, u, u_exact);
      CG1_FORMS::Form_norm2 l2(mesh, u, u_exact);
      double l1_norm = dolfin::assemble(l1);
      double l2_norm = dolfin::assemble(l2);
    
      std::cout << std::setprecision(16) << mesh.hmin() << " " <<
                                            l1_norm << " " <<
                                            l2_norm << " " <<
                                            num_iters << std::endl;

      file << std::setprecision(16) << mesh.hmin() << " " <<
                                       l1_norm << " " <<
                                       l2_norm << " " <<
                                       num_iters << std::endl;
      file.close();
    }
  }
  //----------------------------------------------------------------------------
}

#endif // _TEST_H_

