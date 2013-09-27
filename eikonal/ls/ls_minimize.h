#ifndef _LS_MINIMIZE_H_
#define _LS_MINIMIZE_H_

/*
  Local solvers based on minimizing distance + travel time functionals.
*/

#include <vector>
#include <gsl/gsl_min.h>

namespace eikonal
{
  class LsData
  {
  public:
    LsData(const std::vector<double>& _u_point, const double _u_value, 
           const std::vector<double>& _k_points, 
           const std::vector<double>& _k_values) : 
           u_point(_u_point), u_value(_u_value), k_points(_k_points),
           k_values(_k_values) { }

  public:
    const std::vector<double>& u_point;
    const double u_value;
    const std::vector<double>& k_points;
    const std::vector<double>& k_values;
  };

  double ls_lin_2d_f(double t, void* _params);

  class LinMinimizer2d
  {
    public:
      LinMinimizer2d(gsl_function* _F, double a, double b, double tol, 
                     std::size_t _max_iter);
      ~LinMinimizer2d();
      int find_minimizer(double& t) const;
      double eval(double x) const;
    
    private:
      gsl_function* F;
      const gsl_min_fminimizer_type *T;
      gsl_min_fminimizer *s;

      double lower_bound;
      double upper_bound;
      double tolerance;
      std::size_t max_iter;
  };

  // Given triangle ABC with values of u known at A, B compute value of u at C
  // u_point = C, u_value is guess for the solution
  // k_points = [A, B], k_values = [u[A], u[B]]
  double ls_minimize_2d(const std::vector<double>& u_point,
                        const double u_value,
                        const std::vector<double>& k_points,
                        const std::vector<double>& k_values);

    
}

#endif // _LS_MINIMIZE_H_
