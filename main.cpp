#include "eikonal.h"
#include "test/CG1_VECTOR.h"
#include "test/CG1_FORMS.h"
#include <dolfin.h>
#include <cstdlib>

using namespace eikonal;

// all_linear_tests(atoi(argv[1]), atoi(argv[2]));

int main(int argc, char* argv[])
{
  for(std::size_t i = 3; i < 8; i++)
  {
    std::size_t N = (std::size_t)pow(2, i);
    dolfin::UnitSquareMesh mesh(N, N);
    CG1::FunctionSpace V(mesh);
   
    dolfin::Function u(V);
    dolfin::Function du_dx(V);
    dolfin::Function du_dy(V);
    
    dolfin::Function exact_u(V);
    dolfin::Function exact_du_dx(V);
    dolfin::Function exact_du_dy(V);
   
    std::set<dolfin::la_index> fixed_dofs;

    double _P[2] = {0.5, 0.5}; std::vector<double> P(_P, _P+2);
    MyPoint point(P);
    Problem problem(point);
    problem.init(fixed_dofs, u, du_dx, du_dy);
    problem.exact_solution(exact_u, exact_du_dx, exact_du_dy);

    HermiteSolver solver(V);
    std::size_t precision = 2;
    solver.solve(u, du_dx, du_dy, fixed_dofs, precision);

    dolfin::plot(u);
    dolfin::interactive(true);

    CG1_VECTOR::FunctionSpace VV(mesh);
    dolfin::Function du(VV);
    std::vector<double> du_values(VV.dim());

    boost::shared_ptr<dolfin::GenericVector> dx_vector = du_dx.vector();
    boost::shared_ptr<dolfin::GenericVector> dy_vector = du_dy.vector();
    for(std::size_t i = 0; i < du_values.size()/2; i++)
    {
      du_values[2*i] = (*dx_vector)[i];
      du_values[2*i+1] = (*dy_vector)[i];
    }
   
    du.vector()->set_local(du_values);
    //dolfin::plot(du);
    //#dolfin::interactive(true);

    CG1_FORMS::Form_norm1 l1(mesh, u, exact_u);
    CG1_FORMS::Form_norm2 l2(mesh, u, exact_u);
    double l1_norm = dolfin::assemble(l1);
    double l2_norm = dolfin::assemble(l2);
    
    std::cout << mesh.hmin() << " " << l1_norm << " " << l2_norm << std::endl;
  } 
  return 0;
}

