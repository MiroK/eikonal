#include "ls_minimize.h"
#include "la/la_loop.h"
#include "la/la_common.h"
#include <cassert>
#include <boost/math/tools/minima.hpp>
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
    Linear2DFunctor foo(A, B, C, u_A, u_B);

    //for(size_t i = 0; i < 10; i++)
    //  std::cout << foo(i/10.) << std::endl;

    pair<double, double> t_ft = 
    boost::math::tools::brent_find_minima(foo, 0., 1., 32);

    double u_ = t_ft.second;
    return u_ < u_C ? u_ : u_C;
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
  Linear2DFunctor::Linear2DFunctor(const std::vector<double>& _A,
                                   const std::vector<double>& _B,
                                   const std::vector<double>& _C,
                                   const double& _u_A, const double& _u_B) :
                            A(&_A), B(&_B), C(&_C), u_A(&_u_A), u_B(&_u_B) { }
  //---------------------------------------------------------------------------

  double Linear2DFunctor::operator()(double x)
  {
    const vector<double> P = *A*(1-x) + *B*x;
    const double f = *u_A*(1-x) + *u_B*x + norm(*C-P, 2);
    print(P);
    return f;
  }
  //---------------------------------------------------------------------------
}
