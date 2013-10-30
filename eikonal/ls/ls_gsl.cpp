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

  double hermite_f(double t, void *params)
  {
    HermiteData* data = static_cast<HermiteData*>(params);
    const std::vector<double>& A(data->A), B(data->B), C(data->C),
                             grad_u_A(data->grad_u_A), grad_u_B(data->grad_u_B);
    const double u_A(data->u_A), u_B(data->u_B);

    const vector<double> P = A*(1.-t) + B*t;
    const double f =
    u_A*pow(1-t, 3) + 3*(u_A + dot(B-A, grad_u_A)/3.)*t*(1-t)*(1-t)
    + u_B*pow(t, 3.) + 3*(u_B + dot(A-B, grad_u_B)/3.)*t*t*(1-t)
    + norm(C-P, 2);
    return f;
  }
  //---------------------------------------------------------------------------

  double hermite_df(double t, void *params)
  {
    HermiteData* data = static_cast<HermiteData*>(params);
    const std::vector<double>& A(data->A), B(data->B), C(data->C),
                             grad_u_A(data->grad_u_A), grad_u_B(data->grad_u_B);
    const double u_A(data->u_A), u_B(data->u_B);

    const vector<double> P = A*(1.-t) + B*t;
    const double df = 3*u_B*t*t + 3*(u_B + dot(A-B, grad_u_B)/3.)*t*(2-3*t)
             - 3*u_A*pow(1-t, 2) + 3*(u_A + dot(B-A, grad_u_A)/3.)*(1-t)*(1-3*t)
             + dot(C-P, A-B)/norm(C-P, 2);
    return df;
  }
  //---------------------------------------------------------------------------

  double hermite_ddf(double t, void *params)
  {
    HermiteData* data = static_cast<HermiteData*>(params);
    const std::vector<double>& A(data->A), B(data->B), C(data->C),
                             grad_u_A(data->grad_u_A), grad_u_B(data->grad_u_B);
    const double u_A(data->u_A), u_B(data->u_B);

    const vector<double> P = A*(1.-t) + B*t;

    const double nAB = norm(A-B);
    const double nCP = norm(C-P);
    const double _dot = dot(C-P,A-B);

    const double ddf = 6*u_B*t + 3*(u_B + dot(A-B, grad_u_B)/3.)*(2-6*t)
             + 6*u_A*(1-t) + 3*(u_A + dot(B-A, grad_u_A)/3.)*(-4+6*t)
             + (nAB*nAB*nCP*nCP - _dot*_dot)/nCP/nCP/nCP;
    return ddf;
  }
  //---------------------------------------------------------------------------

  void hermite_df_ddf(double t, void* params, double* y, double* dy)
  {
    *y = hermite_df(t, params);
    *dy = hermite_ddf(t, params);
  }
  //---------------------------------------------------------------------------
  
  std::pair<double, double> 
  hermite_solver(const std::vector<double>& A,
                 const std::vector<double>& B,
                 const std::vector<double>& C,
                 const double u_A, const double u_B, const double u_C,
                 const std::vector<double>& grad_u_A,
                 const std::vector<double>& grad_u_B,
                 std::vector<double>& grad_u_C,
                 const std::size_t max_iter, const std::size_t precision,
                 std::size_t& n_calls)
  {
    // call the newton solver to get the initial guess
    std::pair<double, double> t_ft =
    linear_solver(A, B, C, u_A, u_B, u_C, max_iter, precision, n_calls);

    double t_guess = t_ft.first;
    if(t_guess < 0)// linear did poor job and returned -1, compute the guess
    {
      const double a = norm(C - B, 2);
      const double b = norm(C - A, 2);
      t_guess = (u_A+b) <= (u_B+a) ? 0 : 1;
    }

    // set up the hermite problem
    HermiteData data(A, B, C, u_A, u_B, grad_u_A, grad_u_B);
    gsl_function_fdf F;
    F.f = &hermite_df;
    F.df = &hermite_ddf;
    F.fdf = &hermite_df_ddf;
    F.params = &data; 

    std::size_t newton_calls;
    double t = newton_solver(&F, t_guess, max_iter, precision, newton_calls);
    n_calls += newton_calls; // total of linear and hermite
    
    std::cout << "t " << t << std::endl;

    // check the results
    if(t >= 0 and t <= 1)
    {
      // if smaller value set also the gradient
      double u_ = hermite_f(t, &data);
      if(u_ < u_C)
      {
        std::vector<double> I = A*(1-t) + B*t;
        grad_u_C = (C - I)/norm(C-I, 2);
        
        return std::make_pair<double, double>(t, u_);
      }
      else
      {
        // leave the gradient alone
        return std::make_pair<double, double>(t, u_C);
      }
    }
    else
    {
      // leave the gradient alone
      return std::make_pair<double, double>(-1., u_C);
    }
  }
 //---------------------------------------------------------------------------
}
