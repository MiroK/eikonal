#include "ls_minimize.h"
#include "la/la_loop.h"
#include <cassert>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;

namespace eikonal
{
  
  double ls_minimize_2d(const std::vector<double>& u_point,
                        const double u_value,
                        const std::vector<double>& k_points,
                        const std::vector<double>& k_values)
  {
    const size_t n_k_points = 2;
    const size_t dim = 2;
    assert(u_point.size() == dim);
    assert(k_points.size() == n_k_points*dim);
    assert(k_values.size() == n_k_points);

    LsData ls_data(u_point, u_value, k_points, k_values);
    gsl_function F;
    F.function = &ls_lin_2d_f;
    F.params = &ls_data;

    LinMinimizer2d lin_min_2d(&F, 0, 1, 1E-8, 100);
    
    double t;
    int status = lin_min_2d.find_minimizer(t);
    
    return std::min(u_value, lin_min_2d.eval(t)); // include side compare!
    
  }
  //---------------------------------------------------------------------------

  double ls_lin_2d_f(double t, void* _params)
  {
    std::size_t dim = 2;
    LsData* params = static_cast<LsData*>(_params);
    const double u_A = params->k_values[0],
                 u_B = params->k_values[1],
                 u_C = params->u_value;
    
    const double* _A = &((params->k_points)[0*dim]);
    const double* _B = &((params->k_points)[1*dim]);

    const std::vector<double> A(_A, _A + dim), B(_B, _B + dim),
                              C(params->u_point);
    
    const vector<double> P = A*(1-t) + B*t;
    const double f = u_A*(1-t) + u_B*t + norm(C-P, 2);
    return f;  
  }
  //---------------------------------------------------------------------------

  LinMinimizer2d::LinMinimizer2d(gsl_function* _F, double a, double b,
                               double tol, std::size_t _max_iter) : 
    F(_F), lower_bound(a), upper_bound(b), tolerance(tol), max_iter(_max_iter),
    T(gsl_min_fminimizer_brent)
  {
    s = gsl_min_fminimizer_alloc(T);
  }
  //---------------------------------------------------------------------------

  int LinMinimizer2d::find_minimizer(double& t) const
  {
    // set initial guess
    double a = lower_bound, b = upper_bound;
    t = 0.5;
    gsl_min_fminimizer_set(s, F, t, a, b);
    size_t iter = 0;
    int status;
    do
    {
      iter++;
      status = gsl_min_fminimizer_iterate(s);
      t = gsl_min_fminimizer_x_minimum(s);
      a = gsl_min_fminimizer_x_lower(s);
      b = gsl_min_fminimizer_x_upper(s);
      status = gsl_min_test_interval(a, b, tolerance, tolerance);
    }
    while(status == GSL_CONTINUE and iter < max_iter);

    return status;
  }
  //---------------------------------------------------------------------------
      
  double LinMinimizer2d::eval(double x) const
  {
   return GSL_FN_EVAL(F, x);
  }
  //---------------------------------------------------------------------------
  
  LinMinimizer2d::~LinMinimizer2d()
  {
    gsl_min_fminimizer_free(s);
  }
  //---------------------------------------------------------------------------

}
