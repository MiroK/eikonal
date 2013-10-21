#include "HermiteTest.h"
#include "Problem.h"
#include "MeshGenerator.h"
#include "test_common.h"
#include "CG1.h"
#include "CG1_FORMS.h"
#include "CG1_VECTOR.h"
#include "gs/HermiteSolver.h"
#include <dolfin/mesh/Mesh.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/io/File.h>
#include <dolfin/plot/plot.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>

using namespace dolfin;

namespace eikonal
{

  int hermite_test(const Problem& problem, MeshGenerator& mesh_gen,
                   std::size_t precision, bool plot_on)
  {
    // get the file names
    std::ostringstream help;
    help << problem.name().c_str() << "_" << HermiteSolver::name.c_str();
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
                 "   time     | obtuse/cells |   n_iters    |   min_calls   |"
                 "  max_calls |" << std::endl; 
      
    for( ; !mesh_gen.end(); ++mesh_gen)
    {
      // construct a mesh
      boost::shared_ptr<const dolfin::Mesh> mesh = *mesh_gen;
      const std::size_t num_cells = mesh->num_cells();
      const std::size_t n_obtuse_cells = asset_obtuse_cells(*mesh).size();

      // get the file name for solution, exact_solution and error
      help.str("");
      help << "u_" << data_file_name << "_" << num_cells << ".pvd"; 
      std::string u_file_name = help.str(); // u_problem_solver_N.pvd

      help.str("");
      help << "exact_" << data_file_name << "_" << num_cells << ".pvd"; 
      std::string exact_file_name = help.str(); // exact_problem_solver_N.pvd

      help.str("");
      help << "error_" << data_file_name << "_" << num_cells << ".pvd";
      std::string error_file_name = help.str(); // error_problem...

      // get the solution
      std::size_t num_iters, min_calls = 0, max_calls = 0;
      double l1_norm, l2_norm, coo_norm, time;
      
      int status =
      hermite_test(problem, *mesh, precision, num_iters, min_calls,
                   max_calls, l1_norm, l2_norm, coo_norm, time,
                   u_file_name, exact_file_name, error_file_name, plot_on); 

      std::cout.precision(8);
      std::cout.width(14);
      if(max_calls != 1E6 or min_calls != 0) // write for iterative solver
      {
      // write to screen
      std::cout << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " <<
        n_obtuse_cells << "/" << num_cells << " " << num_iters << " " <<
        min_calls << " " << max_calls << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      
      data_file << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " << n_obtuse_cells <<
        " " << num_cells << " " << num_iters << " " << min_calls << 
        " " << max_calls << std::endl;
      
      data_file.close();
      }

      if(max_calls == 1E6 and min_calls == 0) 
      {
      // write to screen
      std::cout << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " <<
        n_obtuse_cells << "/" << num_cells << " " << num_iters << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      
      data_file << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " << n_obtuse_cells <<
        " " << num_cells << " " << num_iters << std::endl;
      
      data_file.close();
      }
    }
    std::cout << std::endl;

  }
  //---------------------------------------------------------------------------

  int hermite_test(const Problem& problem, const dolfin::Mesh& mesh,
                      std::size_t precision, std::size_t& num_iters,
                      std::size_t& min_calls, std::size_t& max_calls,
                      double& l1_norm, double& l2_norm,
                      double& coo_norm, double& time,
                      std::string u_file_name, std::string exact_file_name,
                      std::string error_file_name, bool plot_on)
  {
    CG1::FunctionSpace V(mesh);
    
    // solver variables, result eval
    dolfin::Function u(V);
    dolfin::Function du_dx(V);
    dolfin::Function du_dy(V);
    std::set<dolfin::la_index> fixed_dofs;
    
    dolfin::Function exact_u(V);
    dolfin::Function exact_du_dx(V);
    dolfin::Function exact_du_dy(V);
  
    // set up vars with solver
    problem.init(fixed_dofs, u, du_dx, du_dy);
    problem.exact_solution(exact_u, exact_du_dx, exact_du_dy);
  
    // set up solver and solve
    HermiteSolver solver(V);

    clock_t start = clock();
    solver.solve(u, du_dx, du_dy, fixed_dofs, precision);
    time = (double)(clock() - start)/CLOCKS_PER_SEC;
  
    // get the error in norms
    CG1_FORMS::Form_norm1 l1(mesh, u, exact_u);
    CG1_FORMS::Form_norm2 l2(mesh, u, exact_u);
    l1_norm = dolfin::assemble(l1);
    l2_norm = dolfin::assemble(l2);
    
    // save exact solution now, because later it is used to hold the error
    dolfin::File exact_file(exact_file_name);
    exact_file << exact_u;
    
    *exact_u.vector() -= *u.vector();
    exact_u.vector()->abs();
    coo_norm = exact_u.vector()->max();

    // if this is an iterative solver extract min/max number of calls to eval
    min_calls = solver.min_calls;
    max_calls = solver.max_calls;

    // save solution
    dolfin::File u_file(u_file_name);
    u_file << u;

    // save error
    dolfin::File error_file(error_file_name);
    error_file << exact_u;

    // plot
    if(plot_on)
    {
      // get the vector ready
      CG1_VECTOR::FunctionSpace VV(mesh);
      
      dolfin::Function du(VV);
      dolfin::Function exact_du(VV);
      std::vector<double> du_values(VV.dim());
      std::vector<double> exact_du_values(VV.dim());
      
      boost::shared_ptr<dolfin::GenericVector> du_dx_vector = du_dx.vector();
      boost::shared_ptr<dolfin::GenericVector> du_dy_vector = du_dy.vector();
      boost::shared_ptr<dolfin::GenericVector>
      exact_du_dx_vector = exact_du_dx.vector();
      boost::shared_ptr<dolfin::GenericVector>
      exact_du_dy_vector = exact_du_dy.vector();
      
      for(std::size_t i = 0; i < du_values.size()/2; i++)
      {
        du_values[2*i] = (*du_dx_vector)[i];
        du_values[2*i+1] = (*du_dy_vector)[i];
        exact_du_values[2*i] = (*exact_du_dx_vector)[i];
        exact_du_values[2*i+1] = (*exact_du_dy_vector)[i];
      }
      du.vector()->set_local(du_values);
      exact_du.vector()->set_local(exact_du_values);

      dolfin::plot(u);
      dolfin::plot(exact_u); // this is the error
      dolfin::interactive(true);
    }
  }
  //---------------------------------------------------------------------------
}
