#ifndef _LS_GSL_H_
#define _LS_GSL_H_

#include <gsl/gsl_math.h>
#include <vector>
#include <utility>

namespace eikonal
{
  // find the root of F, from init guess t_guess,
  // precision is std::numeric_limits<double>::digits10/precision
  // 1 gives 10^-{15}, n_calls holds number of iterations used to find minima
  double
  newton_solver(gsl_function_fdf* F,
                const double t_guess,
                const std::size_t max_iter,
                const std::size_t precision,
                std::size_t& n_calls);

  // LINEAR INTERPOLATION
  class LinearData
  {
  // wrapper of data used by linear solver
  public:
    LinearData(const std::vector<double>& _A, const std::vector<double>& _B,
               const std::vector<double>& _C, const double _u_A,
               const double _u_B) : A(_A), B(_B), C(_C), u_A(_u_A), u_B(_u_B) { }

  public:
    const std::vector<double>& A;
    const std::vector<double>& B;
    const std::vector<double>& C;
    const double u_A;
    const double u_B;
  };

  // evaluate the linearly interpolated distance function
  double linear_f(double t, void *params);

  // evaluate the derivative of the linearly interpolated distance function
  double linear_df(double t, void *params);

  // evaluate the second derivative of the linearly interpolated distance function
  double linear_ddf(double t, void *params);

  // combine df and ddf
  void linear_df_ddf(double t, void* params, double *y, double *dy);

  // given input arguments, construct the interpolant and find its minimal value
  // returns t* where minima obtained and f(t*), n_calls holds number of
  // iterations used in newton solver
  std::pair<double, double>
  linear_solver(const std::vector<double>& A, const std::vector<double>& B,
                const std::vector<double>& C, const double u_A,
                const double u_B, const double u_C, const std::size_t max_iter,
                const std::size_t precision, std::size_t& n_calls);
  
  // HERMITE INTERPOLATION
  class HermiteData : public LinearData
  {
    // wrapper of data used by Hermite solver
  public:
    HermiteData(const std::vector<double>& _A, const std::vector<double>& _B,
                const std::vector<double>& _C, const double _u_A,
                const double _u_B, const std::vector<double>& _grad_u_A,
                const std::vector<double>& _grad_u_B) : 
            LinearData(_A, _B, _C, _u_A, _u_B), grad_u_A(_grad_u_A),
            grad_u_B(_grad_u_B){ }
  public:  
    const std::vector<double>& grad_u_A;
    const std::vector<double>& grad_u_B;

  };

  // evaluate the distance function interpolated with Hermite cubic polynomial 
  double hermite_f(double t, void *params);

  // evaluate the derivative of the distance function
  double hermite_df(double t, void *params);

  // evaluate the second derivative
  double hermite_ddf(double t, void *params);

  // evaluate second and first derivatives
  void hermite_df_ddf(double t, void* params, double* y, double* dy);
  
  // construct the hermite interpolant from data and find its minima
  // returs t*, f(t*), n_calls counts calls made to newton solver, grad_u_C
  // hold approximation of the gradient
  std::pair<double, double> 
  hermite_solver(const std::vector<double>& A,
                 const std::vector<double>& B,
                 const std::vector<double>& C,
                 const double u_A, const double u_B, const double u_C,
                 const std::vector<double>& grad_u_A,
                 const std::vector<double>& grad_u_B,
                 std::vector<double>& grad_u_C,
                 const std::size_t max_iter, const std::size_t precision,
                 std::size_t& n_calls);
}

#endif // _LS_GSL_H_
