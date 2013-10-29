#include "ls_gsl.h"
#include "la/la_loop.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <limits>

#include <iostream>

using namespace std;

namespace eikonal
{
  double
  newton_solver(gsl_function_fdf* F,
                const double t_guess,
                const std::size_t max_iter,
                const std::size_t precision,
                std::size_t& n_calls)
  {
    int status;
    int iter = 0;
    double tol = pow(10., -1.*std::numeric_limits<double>::digits10/precision);
    
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;

    double t = t_guess, t0 = t_guess;
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, F, t);

    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      t0 = t;
      t = gsl_root_fdfsolver_root (s);
      double res = GSL_FN_FDF_EVAL_F(F, t);
      status = gsl_root_test_residual (res, tol);
    } while(status == GSL_CONTINUE and iter < max_iter);
    
    gsl_root_fdfsolver_free (s);
    
    n_calls = iter;
    return t;
  }
  //-----------------------------------------------------------------------------

  double linear_f(double t, void *params)
  {
    LinearData* data = static_cast<LinearData*>(params);
    const vector<double> P = (data->A)*(1.-t) + (data->B)*t;
    const double f = (data->u_A)*(1.-t) + (data->u_B)*t + norm((data->C)-P, 2);
    return f;
  }
  //----------------------------------------------------------------------------

  double linear_df(double t, void *params)
  {
    LinearData* data = static_cast<LinearData*>(params);
    const vector<double> P = (data->A)*(1.-t) + (data->B)*t;
    const double df = -(data->u_A) + (data->u_B) +
                     dot((data->C)-P, (data->A)-(data->B))/norm((data->C)-P, 2);
    return df;
  }
  //----------------------------------------------------------------------------

  double linear_ddf(double t, void *params)
  {
    LinearData* data = static_cast<LinearData*>(params);
    const vector<double> P = (data->A)*(1.-t) + (data->B)*t;

    const double nAB = norm((data->A)-(data->B), 2);
    const double nCP = norm((data->C)-P, 2);
    const double _dot = dot((data->C)-P,(data->A)-(data->B));

    const double ddf = (nAB*nAB*nCP*nCP - _dot*_dot)/nCP/nCP/nCP;
    return ddf;
  }
  //----------------------------------------------------------------------------

  void linear_df_ddf(double t, void* params, double *y, double *dy)
  {
    *y = linear_df(t, params);
    *dy = linear_ddf(t, params);
  }
  //---------------------------------------------------------------------------

  std::pair<double, double> 
  linear_solver(const std::vector<double>& A, const std::vector<double>& B,
                const std::vector<double>& C, const double u_A,
                const double u_B, const double u_C, const std::size_t max_iter,
                const std::size_t precision, std::size_t& n_calls)
  {
    // wrap the data and setup the F
    LinearData data(A, B, C, u_A, u_B);
    gsl_function_fdf F;
    F.f = &linear_df;
    F.df = &linear_ddf;
    F.fdf = &linear_df_ddf;
    F.params = &data; 

    // make an initial guess
    const double a = norm(C - B, 2);
    const double b = norm(C - A, 2);
    double t_guess = (u_A+b) <= (u_B+a) ? 0 : 1;

    double t = newton_solver(&F, t_guess, max_iter, precision, n_calls);
    // see about t
    if(t >= 0 and t <= 1)
    {
      double ft = linear_f(t, &data);
      return std::make_pair<double, double>(t, ft < u_C ? ft : u_C);
    }
    else
    {
      return std::make_pair<double, double>(-1, u_C);
    }
  }
  //----------------------------------------------------------------------------
}
