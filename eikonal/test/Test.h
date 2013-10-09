#ifndef _TEST_H_
#define _TEST_H_

#include "Problem.h"
#include "gs/Solver.h"
#include "CG1.h"
#include "CG1_FORMS.h"
#include <dolfin/generation/UnitSquareMesh.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/io/File.h>
#include <dolfin/plot/plot.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

/*
  Perform convergence test of the Eikonal solver on different problems.
*/

namespace eikonal
{
  // use CG1 in 2D to test Eikonal solvers on meshes of type UnitSquareMesh
  // with 2^i_min, 2^i_max elements, enable ploting every solution to screen
  template<typename T> int linear_2D_test(const Problem& problem,
                                          const std::size_t i_min=3,
                                          const std::size_t i_max=5,
                                          bool plot_on=false);

  // use CG1 in 2D to test Eikonal solvers on meshes loaded from files
  // enable ploting every solution to screen
  template<typename T> int linear_2D_test(const Problem& problem,
                                  const std::vector<std::string>& mesh_names,
                                  bool plot_on=false);

  // use CG1 in 2D to test Eikonal solver on a given mesh
  template<typename T> int linear_2D_test(const Problem& problem,
                                          const dolfin::Mesh& mesh,
                                          std::size_t& num_iters,
                                          double& l1_norm,
                                          double& l2_norm,
                                          std::string u_file_name,
                                          bool plot_on=false);
}

// implementation

namespace eikonal
{

  template<typename T> int linear_2D_test(const Problem& problem,
                                          const std::size_t i_min,
                                          const std::size_t i_max,
                                          bool plot_on=false)
  {
    // get the file names
    std::ostringstream help;
    help << problem.name().c_str() << "_" << T::name.c_str();
    std::string file_name_base = help.str();       // problem_solver

    std::string mesh_type("f"); // uses Fenics generated meshes
    help << "_" << mesh_type.c_str() << ".dat";
    std::string data_file_name = help.str();    // problem_solver_f.dat

    std::ofstream data_file;
    data_file.open(data_file_name.c_str(), std::ios::out);
    data_file.close();
    for(std::size_t i = i_min; i < i_max; i++)
    {
      // construct a mesh
      const std::size_t N = (std::size_t)pow(2, i);
      dolfin::UnitSquareMesh mesh(N, N);
      const std::size_t num_cells = mesh.num_cells();

      // get the file name for solution
      help.str("");
      help << "u_" << file_name_base << "_" << num_cells << ".pvd"; 
      std::string u_file_name = help.str(); // u_problem_solver_N.pvd

      // get the solution
      std::size_t num_iters;
      double l1_norm, l2_norm;
      
      int status =
      linear_2D_test<T>(problem, mesh, num_iters, l1_norm, l2_norm,
                        u_file_name, plot_on); 

      // write to screen
      std::cout << std::setprecision(16) << mesh.hmin() << " " <<
                                            l1_norm << " " <<
                                            l2_norm << " " <<
                                            num_iters << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      data_file << std::setprecision(16) << mesh.hmin() << " " <<
                                            l1_norm << " " <<
                                            l2_norm << " " <<
                                            num_iters << std::endl;
      data_file.close();
    }
  }
  //----------------------------------------------------------------------------
  
  template<typename T> int linear_2D_test(const Problem& problem,
                                  const std::vector<std::string>& mesh_names,
                                  bool plot_on=false)
  {
    // get the file names
    std::ostringstream help;
    help << problem.name().c_str() << "_" << T::name.c_str();
    std::string file_name_base = help.str();       // problem_solver

    std::string mesh_type("g"); // uses Gmsh generated meshes
    help << "_" << mesh_type.c_str() << ".dat";
    std::string data_file_name = help.str();    // problem_solver_f.dat

    std::ofstream data_file;
    data_file.open(data_file_name.c_str(), std::ios::out);
    data_file.close();

    std::vector<std::string>::const_iterator mesh_name = mesh_names.begin();
    for( ; mesh_name != mesh_names.end(); mesh_name++)
    {
      // construct a mesh
      dolfin::UnitSquareMesh mesh(*mesh_name);
      std::size_t num_cells = mesh.num_cells();

      // get the file name for solution
      help.str("");
      help << "u_" << file_name_base << "_" << num_cells << ".pvd"; 
      std::string u_file_name = help.str(); // u_problem_solver_N.pvd

      // get the solution
      std::size_t num_iters;
      double l1_norm, l2_norm;
      
      int status =
      linear_2D_test<T>(problem, mesh, num_iters, l1_norm, l2_norm,
                        u_file_name, plot_on); 

      // write to screen
      std::cout << std::setprecision(16) << mesh.hmin() << " " <<
                                            l1_norm << " " <<
                                            l2_norm << " " <<
                                            num_iters << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      data_file << std::setprecision(16) << mesh.hmin() << " " <<
                                            l1_norm << " " <<
                                            l2_norm << " " <<
                                            num_iters << std::endl;
      data_file.close();
    }

  }
  //----------------------------------------------------------------------------

  template<typename T> int linear_2D_test(const Problem& problem,
                                          const dolfin::Mesh& mesh,
                                          std::size_t& num_iters,
                                          double& l1_norm,
                                          double& l2_norm,
                                          std::string u_file_name,
                                          bool plot_on=false)
  {
      CG1::FunctionSpace V(mesh);

      dolfin::Function u(V);
      dolfin::Function u_exact(V);
      std::set<dolfin::la_index> fixed_dofs;
      
      // set boundary condition, u, and get the exact solution
      problem.init(fixed_dofs, u);
      problem.exact_solution(u_exact);

      // create solver and compute the solution
      T solver(V);
      num_iters = solver.solve(u, fixed_dofs); 

      // get the error of in L1, L2 norms
      CG1_FORMS::Form_norm1 l1(mesh, u, u_exact);
      CG1_FORMS::Form_norm2 l2(mesh, u, u_exact);
      l1_norm = dolfin::assemble(l1);
      l2_norm = dolfin::assemble(l2);

      // save solution
      dolfin::File u_file(u_file_name);
      u_file << u;

      // plot
      if(plot_on)
      {
        dolfin::plot(u);
        dolfin::plot(u_exact);
        dolfin::interactive(true);
      }
  }
  //----------------------------------------------------------------------------
}

#endif // _TEST_H_

