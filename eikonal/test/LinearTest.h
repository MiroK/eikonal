#ifndef _LINEAR_TEST_H_
#define _LINEAR_TEST_H_

#include "Problem.h"
#include "MeshGenerator.h"
#include "CG1.h"
#include "CG1_FORMS.h"
#include "CG1_VECTOR.h"
#include "CG1_VECTOR_FORMS.h"
#include "test_common.h"
#include "gs/Solver.h"
#include <dolfin/fem/assemble.h>
#include <dolfin/io/File.h>
#include <dolfin/plot/plot.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubMesh.h>
#include <dolfin/function/Constant.h>
#include <dolfin/fem/solve.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/solve.h>
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
  // for iterative solver define convergence tolerance as 
  // std::num_limits<double>::digits/precision
  template<typename T> int linear_2D_test(Problem& problem,
                                          MeshGenerator& mesh_gen,
                                          std::size_t precision,
                                          bool plot_on=false);

  // use CG1 in 2D to test Eikonal solver on a given mesh
  // return number of iterations of the sweeping, L^1,L^2 and C^{oo} norms of the
  // error + thec comp time and save the solution to file under u_file_name
  template<typename T> int linear_2D_test(Problem& problem,
                                          const dolfin::Mesh& mesh,
                                          std::size_t precision,
                                          std::size_t& num_iters,
                                          std::size_t& min_calls,
                                          std::size_t& max_calls,
                                          double& l1_norm,
                                          double& l2_norm,
                                          double& band_l1_norm,
                                          double& band_l2_norm,
                                          double& v_l2_norm,
                                          double& band_v_l2_norm,
                                          double& coo_norm,
                                          double& time,
                                          std::string u_file_name,
                                          std::string exact_file_name,
                                          std::string error_file_name,
                                          bool plot_on=false);
}

// implementation

namespace eikonal
{
  template<typename T> int linear_2D_test(Problem& problem,
                                          MeshGenerator& mesh_gen,
                                          std::size_t precision,
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
      double l1_norm, l2_norm, band_l1_norm, band_l2_norm,
             v_l2_norm, band_v_l2_norm, coo_norm, L1_norm, time;
      
      int status =
      linear_2D_test<T>(problem, *mesh, precision, num_iters, min_calls, max_calls,
                        l1_norm, l2_norm, band_l1_norm, band_l2_norm,
                        v_l2_norm, band_v_l2_norm, coo_norm, L1_norm,
                        time, u_file_name, exact_file_name, error_file_name, plot_on); 

      std::cout.precision(8);
      std::cout.width(14);
      if(max_calls != 1E6 or min_calls != 0) // write for iterative solver
      {
      // write to screen
      std::cout << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " <<
        n_obtuse_cells << "/" << num_cells << " " << num_iters << " " <<
        min_calls << " " << max_calls << " " << band_l1_norm << " "
        << band_l2_norm << " " << v_l2_norm << " " << band_v_l2_norm << " " 
        << L1_norm << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      
      data_file << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " << n_obtuse_cells <<
        " " << num_cells << " " << num_iters << " " << min_calls << 
        " " << max_calls << " " << band_l1_norm << " " << band_l2_norm <<
        " " << v_l2_norm << " " << band_v_l2_norm << " " << L1_norm << std::endl;
      
      data_file.close();
      }

      if(max_calls == 1E6 and min_calls == 0) 
      {
      // write to screen
      std::cout << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " <<
        n_obtuse_cells << "/" << num_cells << " " << num_iters <<
        " " << band_l1_norm << " " << band_l2_norm <<
        " " << v_l2_norm << " " << band_v_l2_norm << " " << L1_norm << std::endl;
      
      // write to text file
      data_file.open(data_file_name.c_str(), std::ios::app);
      
      data_file << std::scientific << mesh->hmin() << " " << l1_norm << " " <<
        l2_norm << " " << coo_norm << " " << time << " " << n_obtuse_cells <<
        " " << num_cells << " " << num_iters << 
        " " << band_l1_norm << " " << band_l2_norm <<
        " " << v_l2_norm << " " << band_v_l2_norm << " " << L1_norm << std::endl;
      
      data_file.close();
      }
    }
    std::cout << std::endl;
  }
  //----------------------------------------------------------------------------
  
  template<typename T> int linear_2D_test(Problem& problem,
                                          const dolfin::Mesh& mesh,
                                          std::size_t precision, 
                                          std::size_t& num_iters,
                                          std::size_t& min_calls,
                                          std::size_t& max_calls,
                                          double& l1_norm,
                                          double& l2_norm,
                                          double& band_l1_norm,
                                          double& band_l2_norm,
                                          double& v_l2_norm,
                                          double& band_v_l2_norm,
                                          double& coo_norm,
                                          double& L1_norm,
                                          double& time,
                                          std::string u_file_name,
                                          std::string exact_file_name,
                                          std::string error_file_name,
                                          bool plot_on=false)
  {
      CG1::FunctionSpace V(mesh);

      dolfin::Function u(V);
      dolfin::Function exact_du_dx(V); // components of exact
      dolfin::Function exact_du_dy(V);
      dolfin::Function u_exact(V);
      std::set<dolfin::la_index> fixed_dofs;
      
      // set boundary condition, u, and get the exact solution and band
      problem.init(fixed_dofs, u);
      problem.exact_solution(u_exact, exact_du_dx, exact_du_dy);
      dolfin::MeshFunction<std::size_t> band = problem.get_band(u, 3);

      // create solver and compute the sol //TODO do via call!ution
      T solver(V);
      // time the sweeping
      clock_t start = clock();
      num_iters = solver.solve(u, fixed_dofs, precision); 
      time = (double)(clock() - start)/CLOCKS_PER_SEC;

      // get the error in norms
      CG1_FORMS::Form_norm1 l1(mesh, u, u_exact);
      CG1_FORMS::Form_norm2 l2(mesh, u, u_exact);
      l1_norm = dolfin::assemble(l1);
      l2_norm = sqrt(dolfin::assemble(l2));

      // get the area of the band and norms1 on the submesh
      dolfin::SubMesh band_mesh(mesh, band, 1);
      CG1_FORMS::Form_norm1 l1_band(band_mesh, u, u_exact);
      CG1_FORMS::Form_norm2 l2_band(band_mesh, u, u_exact);
      dolfin::Constant one(1.);
      CG1_FORMS::Form_area area_band(band_mesh, one);

      double area = dolfin::assemble(area_band);
      band_l1_norm = dolfin::assemble(l1_band)/area;
      band_l2_norm = sqrt(dolfin::assemble(l2_band)/area);
     
      // get the vector
      CG1_VECTOR::FunctionSpace VV(mesh);
      
      dolfin::Function du(VV), exact_du(VV);
      std::vector<double> exact_du_values(VV.dim());
      
      boost::shared_ptr<dolfin::GenericVector>
      exact_du_dx_vector = exact_du_dx.vector(),
      exact_du_dy_vector = exact_du_dy.vector();
    
      // set the values of exact
      for(std::size_t i = 0; i < exact_du_values.size()/2; i++)
      {
        exact_du_values[2*i] = (*exact_du_dx_vector)[i];
        exact_du_values[2*i+1] = (*exact_du_dy_vector)[i];
      }
      exact_du.vector()->set_local(exact_du_values);

      // get gradient of numerical solution
      CG1_VECTOR::BilinearForm a(VV, VV);
      CG1_VECTOR::LinearForm L(VV);
      L.f = u;

      dolfin::Matrix A;
      dolfin::Vector b;
      dolfin::assemble_system(A, b, a, L);
      dolfin::solve(A, *du.vector(), b, "cg");

      CG1_VECTOR_FORMS::Form_norm vector_l2_norm(mesh, du, exact_du);
      v_l2_norm = sqrt(dolfin::assemble(vector_l2_norm));
      
      // get the area of the band and norms1 on the submesh
      CG1_VECTOR_FORMS::Form_norm vector_l2_band_norm(band_mesh, du, exact_du);
      band_v_l2_norm = sqrt(dolfin::assemble(vector_l2_band_norm)/area);

      // save exact solution now, because later it is used to hold the error
      dolfin::File exact_file(exact_file_name);
      exact_file << u_exact;
      
      // error vector
      *u_exact.vector() -= *u.vector();
      u_exact.vector()->abs();
      L1_norm = u_exact.vector()->sum()/u_exact.vector()->size();
      coo_norm = u_exact.vector()->max();

      // if this is an iterative solver extract min/max number of calls to eval
      min_calls = solver.min_calls;
      max_calls = solver.max_calls;

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

#endif // _LINEAR_TEST_H_

