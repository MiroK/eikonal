#include "ls_minimize.h"
#include "ls_gsl.h"
#include "la/la_loop.h"
#include "la/la_common.h"
#include <cassert>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <cmath>

using namespace std;

namespace eikonal
{
  std::size_t Linear2DFunctor::n_calls = 0;

  std::pair<double, double>
  linear_brent_2d(const std::vector<double>& A,
                         const std::vector<double>& B,
                         const std::vector<double>& C,
                         const double u_A, const double u_B,
                         const double u_C, std::size_t& n_calls,
                         const std::size_t precision)
  {
    // use boost to minimize via functor
    LinearBrent2D foo(A, B, C, u_A, u_B);
    foo.n_calls = 0;

    int digits = std::numeric_limits<double>::digits/precision;
    pair<double, double> t_ft = 
    boost::math::tools::brent_find_minima(foo, 0., 1., digits);
    
    n_calls = foo.n_calls;

    double u_ = t_ft.second;
    return std::make_pair<double, double>(t_ft.first, u_ < u_C ? u_ : u_C);
  }
  //---------------------------------------------------------------------------
  
  std::pair<double, double>
  linear_newton_2d(const std::vector<double>& A,
                   const std::vector<double>& B,
                   const std::vector<double>& C,
                   const double u_A, const double u_B,
                   const double u_C, std::size_t& n_calls,
                   std::size_t precision)
  {
    // use boost to minimize via functor
    LinearNewton2D foo(A, B, C, u_A, u_B);
    foo.n_calls = 0;
    int digits = std::numeric_limits<double>::digits/precision; 
    double t_min = 0.;
    double t_max = 1.;
    const double a = norm(C - B, 2);
    const double b = norm(C - A, 2);
    double t_guess = (u_A+b) <= (u_B+a) ? 0 : 1;
    double t = boost::math::tools::newton_raphson_iterate(foo,
                                                          t_guess,
                                                          t_min,
                                                          t_max,
                                                          digits);
    n_calls = foo.n_calls;

    if(t >= 0 and t <= 1)
    {
      double u_ = foo.eval(t);
      return std::make_pair<double, double>(t, u_ < u_C ? u_ : u_C);
    }
    else
      return std::make_pair<double, double>(-1., u_C);
  }
  //---------------------------------------------------------------------------

  double linear_brent_2d(const std::vector<double>& u_point,
                         const double u_value,
                         const std::vector<double>& k_points,
                         const std::vector<double>& k_values,
                         std::size_t& n_calls,
                         const std::size_t precision)
  {
    const size_t n_k_points = 2;
    const size_t dim = 2;
    assert(u_point.size() == dim);
    assert(k_points.size() == n_k_points*dim);
    assert(k_values.size() == n_k_points);

    // unpack 
    const std::size_t i = k_values[0] <= k_values[1] ? 0 : 1;
    const std::size_t ip1 = (i+1)%2;
    const double u_A = k_values[i], u_B = k_values[ip1], u_C = u_value;
    
    const std::vector<double> A(&k_points[i*dim], &k_points[i*dim] + dim),
                              B(&k_points[ip1*dim], &k_points[ip1*dim] + dim),
                              C(u_point);
    
    // use the verbose interface, getint 
    // t - parameter of the intersect and the value t, return the value
    return linear_brent_2d(A, B, C, u_A, u_B, u_C, n_calls, precision).second;
  }
  //---------------------------------------------------------------------------

  double linear_newton_2d(const std::vector<double>& u_point,
                          const double u_value,
                          const std::vector<double>& k_points,
                          const std::vector<double>& k_values,
                          std::size_t& n_calls,
                          const std::size_t precision)
  {
    const size_t n_k_points = 2;
    const size_t dim = 2;
    assert(u_point.size() == dim);
    assert(k_points.size() == n_k_points*dim);
    assert(k_values.size() == n_k_points);

    // unpack 
    const std::size_t i = k_values[0] <= k_values[1] ? 0 : 1;
    const std::size_t ip1 = (i+1)%2;
    const double u_A = k_values[i], u_B = k_values[ip1], u_C = u_value;
    
    const std::vector<double> A(&k_points[i*dim], &k_points[i*dim] + dim),
                              B(&k_points[ip1*dim], &k_points[ip1*dim] + dim),
                              C(u_point);
    
    // use the verbose interface
    return linear_newton_2d(A, B, C, u_A, u_B, u_C, n_calls, precision).second;
    //return linear_solver(A, B, C, u_A, u_B, u_C, 40, precision).second;
  }
  //---------------------------------------------------------------------------

  std::pair<double, double>
  hermite_newton_2d(const std::vector<double>& A,
                    const std::vector<double>& B,
                    const std::vector<double>& C,
                    const double u_A, const double u_B, const double u_C,
                    const std::vector<double>& grad_u_A,
                    const std::vector<double>& grad_u_B,
                    std::vector<double>& grad_u_C,
                    std::size_t& n_calls,
                    const std::size_t precision)
  {
    // first call the newton solver
    std::pair<double, double> t_ft =
    linear_newton_2d(A, B, C, u_A, u_B, u_C, n_calls, precision);

    // now call hermite using t_ft.first as guess if possible
    HermiteNewton2D foo(A, B, C, u_A, u_B, grad_u_A, grad_u_B);
    
    const double a = norm(C - B, 2);
    const double b = norm(C - A, 2);
   
    std::cout.precision(16);
    //std::cout << "t from newton " << t_ft.first << std::endl; 

    double t_guess;
    if(t_ft.first < 0) // linear solver did a poor job
    {
      t_guess = (u_A+b) <= (u_B+a) ? 0 : 1;
    }
    else
    {
      t_guess = t_ft.first;
    }
      
    foo.n_calls = 0;
    int digits = std::numeric_limits<double>::digits/precision;
    double t_min = 0.;
    double t_max = 1.;
    boost::uintmax_t max_iter = 300;
    double t = boost::math::tools::newton_raphson_iterate(foo,
                            t_guess, t_min, t_max, digits, max_iter);
    n_calls = foo.n_calls;
    
    if(t >= 0 and t <= 1)
    {
      // if smaller value set also the gradient
      double u_ = foo.eval(t);
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
  //---------------------------------------------------------------------------
  }
}

//-----------------------------------------------------------------------------

namespace eikonal
{
  Linear2DFunctor::Linear2DFunctor(const std::vector<double>& _A,
                                   const std::vector<double>& _B,
                                   const std::vector<double>& _C,
                                   const double _u_A, const double _u_B) :
                  A(_A), B(_B), C(_C), u_A(_u_A), u_B(_u_B){ } 
  //---------------------------------------------------------------------------

  LinearBrent2D::LinearBrent2D(const std::vector<double>& _A,
                               const std::vector<double>& _B,
                               const std::vector<double>& _C,
                               const double _u_A, const double _u_B) :
                      Linear2DFunctor(_A, _B, _C, _u_A, _u_B){ }
  //---------------------------------------------------------------------------

  LinearNewton2D::LinearNewton2D(const std::vector<double>& _A,
                                 const std::vector<double>& _B,
                                 const std::vector<double>& _C,
                                 const double _u_A, const double _u_B) : 
                      Linear2DFunctor(_A, _B, _C, _u_A, _u_B) { }
  //---------------------------------------------------------------------------
  
  double LinearBrent2D::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double f = u_A*(1.-x) + u_B*x + norm(C-P, 2);
    n_calls++;
    return f;
  }
  //---------------------------------------------------------------------------
  
  double LinearNewton2D::eval(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double f = u_A*(1.-x) + u_B*x + norm(C-P, 2);
    return f;
  }
  //---------------------------------------------------------------------------

  boost::math::tuple<double, double> LinearNewton2D::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double df = -u_A + u_B + dot(C-P, A-B)/norm(C-P, 2);       //<---
    
    const double nAB = norm(A-B, 2);
    const double nCP = norm(C-P, 2);
    const double _dot = dot(C-P,A-B);

    const double ddf = (nAB*nAB*nCP*nCP - _dot*_dot)/nCP/nCP/nCP;      //<---
    n_calls++;
    return boost::math::make_tuple(df, ddf);
  }
  //---------------------------------------------------------------------------

  HermiteNewton2D::HermiteNewton2D(const std::vector<double>& _A,
                                   const std::vector<double>& _B,
                                   const std::vector<double>& _C,
                                   const double _u_A, const double _u_B,
                                   const std::vector<double>& _grad_u_A,
                                   const std::vector<double>& _grad_u_B)
  : Linear2DFunctor(_A, _B, _C, _u_A, _u_B), grad_u_A(_grad_u_A),
  grad_u_B(_grad_u_B) { }
  //----------------------------------------------------------------------------

  double HermiteNewton2D::eval(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double f =
    u_A*pow(1-x, 3) + 3*(u_A + dot(B-A, grad_u_A)/3.)*x*(1-x)*(1-x)
    + u_B*pow(x, 3.) + 3*(u_B + dot(A-B, grad_u_B)/3.)*x*x*(1-x)
    + norm(C-P, 2);
    return f;
  }
  //----------------------------------------------------------------------------

  boost::math::tuple<double, double> HermiteNewton2D::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double df = 3*u_B*x*x + 3*(u_B + dot(A-B, grad_u_B)/3.)*x*(2-3*x)
             - 3*u_A*pow(1-x, 2) + 3*(u_A + dot(B-A, grad_u_A)/3.)*(1-x)*(1-3*x)
             + dot(C-P, A-B)/norm(C-P, 2);  //<---

    const double nAB = norm(A-B);
    const double nCP = norm(C-P);
    const double _dot = dot(C-P,A-B);

    const double ddf = 6*u_B*x + 3*(u_B + dot(A-B, grad_u_B)/3.)*(2-6*x)
             + 6*u_A*(1-x) + 3*(u_A + dot(B-A, grad_u_A)/3.)*(-4+6*x)
             + (nAB*nAB*nCP*nCP - _dot*_dot)/nCP/nCP/nCP;      //<---

    n_calls++;
    return boost::math::make_tuple(df, ddf);
  }
  //---------------------------------------------------------------------------

  HermiteBrent2D::HermiteBrent2D(const std::vector<double>& _A,
                                   const std::vector<double>& _B,
                                   const std::vector<double>& _C,
                                   const double _u_A, const double _u_B,
                                   const std::vector<double>& _grad_u_A,
                                   const std::vector<double>& _grad_u_B)
  : Linear2DFunctor(_A, _B, _C, _u_A, _u_B), grad_u_A(_grad_u_A),
  grad_u_B(_grad_u_B) { }
  //----------------------------------------------------------------------------

  double HermiteBrent2D::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;
    const double f =
    u_A*pow(1-x, 3) + 3*(u_A + dot(B-A, grad_u_A)/3.)*x*(1-x)*(1-x)
    + u_B*pow(x, 3.) + 3*(u_B + dot(A-B, grad_u_B)/3.)*x*x*(1-x)
    + norm(C-P, 2);
    n_calls++;
    return f;
  }
  //---------------------------------------------------------------------------

  TEST_FUNCTOR::TEST_FUNCTOR(const std::vector<double>& _A,
                                   const std::vector<double>& _B,
                                   const std::vector<double>& _C,
                                   const double _u_A, const double _u_B,
                                   const std::vector<double>& _grad_u_A,
                                   const std::vector<double>& _grad_u_B)
  : Linear2DFunctor(_A, _B, _C, _u_A, _u_B), grad_u_A(_grad_u_A),
  grad_u_B(_grad_u_B) { }
  //----------------------------------------------------------------------------

  boost::math::tuple<double, double> TEST_FUNCTOR::operator()(double x)
  {
    const vector<double> P = A*(1.-x) + B*x;

    const double f =
    u_A*pow(1-x, 3) + 3*(u_A + dot(B-A, grad_u_A)/3.)*x*(1-x)*(1-x)
    + u_B*pow(x, 3.) + 3*(u_B + dot(A-B, grad_u_B)/3.)*x*x*(1-x)
    + norm(C-P, 2);
    
    const double nAB = norm(A-B);
    const double nCP = norm(C-P);
    const double _dot = dot(C-P,A-B);

    const double df = 3*u_B*x*x + 3*(u_B + dot(A-B, grad_u_B)/3.)*x*(2-3*x)
             - 3*u_A*pow(1-x, 2) + 3*(u_A + dot(B-A, grad_u_A)/3.)*(1-x)*(1-3*x)
             + dot(C-P, A-B)/norm(C-P, 2);  //<---

    n_calls++;
    return boost::math::make_tuple(f, df);
  }
  //---------------------------------------------------------------------------
}
