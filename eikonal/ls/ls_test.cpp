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

  double d_g = linear_geometric_2d(A, B, C, u_A, u_B, u_C);
  double e_g = abs(d_g - 1);
  
  double d_m = linear_brent_2d(A, B, C, u_A, u_B, u_C);
  double e_m = abs(d_m - 1);

  double d_n = linear_newton_2d(A, B, C, u_A, u_B, u_C);
  double e_n = abs(d_n - 1);


  std::cout << "C is "; print(C);
  std::cout << "geometric " << d_g << "\t";
  std::cout << "brent " << d_m << "\t";
  std::cout << "newton " << d_n << "\n";
  std::cout << std::setprecision(16) << e_n << "\t"
                                     << e_m << "\t" 
                                     << e_n << std::endl;

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

  double d_g = linear_geometric_2d(A, B, C, u_A, u_B, u_C);
  double e_g = abs(ex - d_g);
  double d_m = linear_brent_2d(A, B, C, u_A, u_B, u_C);
  double e_m = abs(ex - d_m);

  std::cout << "S is "; print(source);
  std::cout << "geometric " << d_g << "\t";
  std::cout << "minimize " << d_m << "\t";
  std::cout << std::setprecision(16) << abs(d_g - d_m) << std::endl;
  std::cout << std::setprecision(16) << e_g << "\t" << e_m << std::endl;
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
