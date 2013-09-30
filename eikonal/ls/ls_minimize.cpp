#include "ls_minimize.h"
#include "la/la_loop.h"
#include "la/la_common.h"
#include <cassert>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <utility>

using namespace std;

namespace eikonal
{
    
  double linear_minimize_2d(const std::vector<double>& A,
                            const std::vector<double>& B,
                            const std::vector<double>& C,
                            const double u_A, const double u_B,
                            const double u_C)
  {
    // use boost to minimize via functor
    LinearBrent2D foo(A, B, C, u_A, u_B);

    //for(size_t i = 0; i < 10; i++)
    //  std::cout << foo(i/10.) << std::endl;

    std::cout << "\t\t\t\t\tf(0), f(1) " << foo(0) << ", " << foo(1) << "\n";
    pair<double, double> t_ft = 
    boost::math::tools::brent_find_minima(foo, 0., 1., 32);

    double u_ = t_ft.second;
    return u_ < u_C ? u_ : u_C;
  }
  //---------------------------------------------------------------------------
  
  double linear_newton_2d(const std::vector<double>& A,
                           const std::vector<double>& B,
                           const std::vector<double>& C,
                           const double u_A, const double u_B,
                           const double u_C)
  {
    // use boost to minimize via functor
    LinearNewton2D foo(A, B, C, u_A, u_B);
    int digits = std::numeric_limits<double>::digits / 4; 
    double t_min = 0.;
    double t_max = 1.;
    double t_guess = 0.5;
    boost::uintmax_t max_iter = 40;
    double t = boost::math::tools::newton_raphson_iterate(foo,
                                                          t_guess,
                                                          t_min,
                                                          t_max,
                                                          digits,
                                                  max_iter);
    std::cout << "t " << t << std::endl;
    if(t >= 0 and t <= 1)
    {
      double u_ = boost::fusion::get<0>(foo(t));
      return u_ < u_C ? u_ : u_C;
    }
    else
      return u_C;
  }
  //---------------------------------------------------------------------------

  double linear_minimize_2d(const std::vector<double>& u_point,
                            const double u_value,
                            const std::vector<double>& k_points,
                            const std::vector<double>& k_values)
  {
    const size_t n_k_points = 2;
    const size_t dim = 2;
    assert(u_point.size() == dim);
    assert(k_points.size() == n_k_points*dim);
    assert(k_values.size() == n_k_points);

    // unpack 
    const double u_A = k_values[0], u_B = k_values[1], u_C = u_value;
    
    const std::vector<double> A(&k_points[0*dim], &k_points[0*dim] + dim),
                              B(&k_points[1*dim], &k_points[1*dim] + dim),
                              C(u_point);
    
    // use the verbose interface
    return linear_minimize_2d(A, B, C, u_A, u_B, u_C);
  }
  //---------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------

namespace eikonal
{
  double LinearBrent2D::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double f = u_A*(1.-x) + u_B*x + norm(C-P, 2);
    print(P);
    std::cout << "x " << u_A << " f(x) " << u_B << std::endl;
    return f;
  }
  //---------------------------------------------------------------------------

  boost::math::tuple<double, double> LinearNewton2D::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double f = u_A*(1.-x) + u_B*x + norm(C-P, 2);
    const double df = -u_A + u_B + dot(C-P, A-B)/norm(C-P, 2);
    std::cout << "new x " << x << " f(x) " << f << 
                                  " df(x) " << df << std::endl;
    return boost::math::make_tuple(f, df);
  }
}
