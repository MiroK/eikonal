#ifndef _TEST_H_
#define _TEST_H_

#include "Problem.h"
#include "MeshGenerator.h"
#include "CG1.h"
#include "CG1_FORMS.h"
#include "test_common.h"
#include "gs/Solver.h"
#include <dolfin/generation/UnitSquareMesh.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/io/File.h>
#include <dolfin/plot/plot.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>

/*
  Perform convergence test of the Eikonal solver on different problems.
*/

namespace eikonal
{
  // use CG1 in 2D to test Eikonal solvers on meshes of type UnitSquareMesh
  template<typename T> int linear_2D_test(const Problem& problem,
                                          MeshGenerator& mesh_gen,
                                          bool plot_on=false);

  // use CG1 in 2D to test Eikonal solver on a given mesh
  // return number of iterations of the sweeping, L^1,L^2 and C^{oo} norms of the
  // error + thec comp time and save the solution to file under u_file_name
  template<typename T> int linear_2D_test(const Problem& problem,
                                          const dolfin::Mesh& mesh,
                                          std::size_t& num_iters,
                                          double& l1_norm,
                                          double& l2_norm,
                                          double& coo_norm,
                                          double& time,
                                          std::string u_file_name,
                                          bool plot_on=false);
}

// implementation

namespace eikonal
{
  template<typename T> int linear_2D_test(const Problem& problem,
                                          MeshGenerator& mesh_gen,
                                          bool plot_on=false)
  {
    // get the file names
    std::ostringstream help;
    help << problem.name().c_str() << "_" << T::name.c_str();
    std::string file_name_base = help.str();       // problem_solver


    std::string mesh_type = mesh_gen.type(); // get the type of mesh from generator
    help << "_" << mesh_type.c_str() << ".dat";
    std::string data_file_name = help.str();    // problem_solver_f.dat
    
    // print identifiers to screen
    std::cout << data_file_name << std::endl;
    
    std::ofstream data_file;
    data_file.open(data_file_name.c_str(), std::ios::out);
    data_file.close();

    // write legend
    std::cout << "    h_min    |    l1_norm  |     l2_norm  |  coo_norm     |"
                 "   time     | obtuse/cells |   n_iters    " << std::endl; 
      
    for( ; !mesh_gen.end(); ++mesh_gen)
    {
      // construct a mesh
      boost::shared_ptr<const dolfin::Mesh> mesh = *mesh_gen;
      const std::size_t num_cells = mesh->num_cells();
      const std::size_t n_obtuse_cells = asset_obtuse_cells(*mesh).size();

      // get the file name for solution and error
      help.str("");
      help << "u_" << data_file_name << "_" << num_cells << ".pvd"; 
      std::string u_file_name = help.str(); // u_problem_solver_N.pvd
      help.str("");
      help << "e_" << data_file_name << "_" << num_cells << ".pvd";
      std::string error_file_name = help.str();

      // get the solution
      std::size_t num_iters;
      double l1_norm, l2_norm, coo_norm, time;
      
      int status =
      linear_2D_test<T>(problem, *mesh, num_iters, l1_norm, l2_norm, coo_norm,
                        time, u_file_name, error_file_name, plot_on); 

      // write to screen
      std::cout.precision(8);
      std::cout.width(14);
      std::cout << std::scientific << mesh->hmin() << " " <<
                                           l1_norm << " " <<
                                           l2_norm << " " <<
                                           coo_norm << " " <<
                                           time << " " <<
                                  n_obtuse_cells << "/" << num_cells << " " <<
                                           num_iters << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      data_file << std::scientific << mesh->hmin() << " " <<
                                           l1_norm << " " <<
                                           l2_norm << " " <<
                                           coo_norm << " " <<
                                           time << " " <<
                                           n_obtuse_cells << " " << num_cells << " " << 
                                           num_iters << std::endl;
      data_file.close();
    }
    std::cout << std::endl;
  }
  //----------------------------------------------------------------------------
  
  template<typename T> int linear_2D_test(const Problem& problem,
                                          const dolfin::Mesh& mesh,
                                          std::size_t& num_iters,
                                          double& l1_norm,
                                          double& l2_norm,
                                          double& coo_norm,
                                          double& time,
                                          std::string u_file_name,
                                          std::string error_file_name,
                                          bool plot_on=false)
  {
      CG1::FunctionSpace V(mesh);

      dolfin::Function u(V);
      dolfin::Function u_exact(V);
      std::set<dolfin::la_index> fixed_dofs;
      
      // set boundary condition, u, and get the exact solution
      problem.init(fixed_dofs, u);
      problem.exact_solution(u_exact);

      dolfin::plot(u_exact);
      dolfin::plot(u);
      dolfin::interactive(true);

      // create solver and compute the solution
      T solver(V);
      // time the sweeping
      clock_t start = clock();
      num_iters = solver.solve(u, fixed_dofs); 
      time = (double)(clock() - start)/CLOCKS_PER_SEC;

      // get the error in norms
      CG1_FORMS::Form_norm1 l1(mesh, u, u_exact);
      CG1_FORMS::Form_norm2 l2(mesh, u, u_exact);
      l1_norm = dolfin::assemble(l1);
      l2_norm = dolfin::assemble(l2);
      *u_exact.vector() -= *u.vector();
      u_exact.vector()->abs();
      coo_norm = u_exact.vector()->max();

      // save solution
      dolfin::File u_file(u_file_name);
      u_file << u;

      // save error
      dolfin::File error_file(error_file_name);
      error_file << u_exact;

      // plot
      if(plot_on)
      {
        dolfin::plot(u);
        dolfin::plot(u_exact); // this is the error
        dolfin::interactive(true);
      }
  }
  //----------------------------------------------------------------------------
}

#endif // _TEST_H_

