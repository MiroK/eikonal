#include "eikonal_ls.h"
#include "la/eikonal_la.h"
#include "cg/eikonal_cg.h"
#include <iostream>
#include <iomanip>

using namespace eikonal;

void local_test_C(double* _C)
{
  std::vector<double> C(_C, _C + 2);
  double _A[2] = {0., 0.}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1., 0.}; std::vector<double> B(_B, _B + 2);

  double u_C = 100;
  double u_A = 0;
  double u_B = 0;
  double _grad_u_A[2] = {0, 1.}; //same as grad_u_B
  std::vector<double> grad_u_A(_grad_u_A, _grad_u_A + 2);
  std::vector<double> grad_u_B(grad_u_A); 
  
  std::size_t n_calls;
  double d_g = linear_geometric_2d(A, B, C, u_A, u_B, u_C);
  double e_g = abs(d_g - 1);
  
  double d_m = linear_brent_2d(A, B, C, u_A, u_B, u_C, n_calls, 2).second;
  double e_m = abs(d_m - 1);

  double d_n = linear_newton_2d(A, B, C, u_A, u_B, u_C, n_calls, 2).second;
  double e_n = abs(d_n - 1);

  std::vector<double> grad_u_C;
  double d_h = hermite_newton_2d(A, B, C, u_A, u_B, u_C, grad_u_A, grad_u_B,
                                 grad_u_C, n_calls, 2).second;
  double e_h = abs(d_h - 1);

  std::cout << "C is "; print(C);
  std::cout << "geometric " << d_g << "\t";
  std::cout << "brent " << d_m << "\t";
  std::cout << "newton " << d_n << "\t";
  std::cout << "hermite " << d_h << "\n";
  std::cout << std::setprecision(16) << e_n << "\t"
                                     << e_m << "\t" 
                                     << e_n << "\t"
                                     << e_h << std::endl;

  std::cout << std::endl;
}
//-----------------------------------------------------------------------------

void local_test_S(double* _S)
{
  double _C[2] = {0.5, 1.}; std::vector<double> C(_C, _C + 2);
  double _A[2] = {0., 0.}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1., 0.}; std::vector<double> B(_B, _B + 2);

  double u_C = 100;
  std::vector<double> source(_S, _S + 2);
  double u_A = point_point(A, source);
  double u_B = point_point(B, source);
  double ex = point_point(C, source);
  const std::vector<double> grad_u_A = (A-source)/norm(A-source);
  const std::vector<double> grad_u_B = (B-source)/norm(B-source);
  std::vector<double> grad_u_C;

  std::cout << "grad_A "; print(grad_u_A);
  std::cout << "grad B "; print(grad_u_B);
  
  std::size_t n_calls;
  double d_g = linear_geometric_2d(A, B, C, u_A, u_B, u_C);
  double e_g = abs(ex - d_g);

  double d_m = linear_solver(A, B, C, u_A, u_B, u_C, 40, 2, n_calls).second;
  double e_m = abs(ex - d_m);
  std::cout << "number of calls in lin Newton solver " <<n_calls << std::endl;
  
  double d_h = hermite_solver(A, B, C, u_A, u_B, u_C, grad_u_A, grad_u_B,
                                 grad_u_C, 40, 2, n_calls).second;
  double e_h = abs(d_h - ex);
  std::cout << "number of calls in herm Newton solver " <<n_calls << std::endl;

  std::cout << "S is "; print(source);
  std::cout << "exact " << ex << "\n";
  std::cout << "geometric " << d_g << " " << e_g << "\n";
  std::cout << "minimize " << d_m << " " << e_m << "\n";
  std::cout << "Hermite " << d_h << " " << e_h << std::endl;
  
  std::cout << "grad_C "; print(grad_u_C);
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
