#ifndef _LS_MINIMIZE_H_
#define _LS_MINIMIZE_H_

/*
  Local solvers based on minimizing distance + travel time functionals.
*/

#include <vector>
#include <boost/math/tools/tuple.hpp>

namespace eikonal
{
  // Given triangle ABC with values of u known at A, B compute value of u at C
  // from guess u_C by monimizeing linear interpolant 
  // return to n_calls number of calls to the internal eval
  // testing interface for Brent method
  double linear_brent_2d(const std::vector<double>& A,
                         const std::vector<double>& B,
                         const std::vector<double>& C,
                         const double u_A, const double u_B, const double u_C,
                         std::size_t& n_calls);
  
  // testing interface for Newton method
  double linear_newton_2d(const std::vector<double>& A,
                          const std::vector<double>& B,
                          const std::vector<double>& C,
                          const double u_A, const double u_B, const double u_C,
                          std::size_t& n_calls);
 
  // global solver interface for Brent method
  // u_point = C, u_value is guess for the solution u_C
  // k_points = [A, B], k_values = [u[A], u[B]]
  double linear_brent_2d(const std::vector<double>& u_point,
                         const double u_value,
                         const std::vector<double>& k_points,
                         const std::vector<double>& k_values,
                         std::size_t& n_calls);
 
  // global solver interface to Newton method
  double linear_newton_2d(const std::vector<double>& u_point,
                          const double u_value,
                          const std::vector<double>& k_points,
                          const std::vector<double>& k_values,
                          std::size_t& n_calls);
}

namespace eikonal
{
  class Linear2DFunctor
  { // parent for linear minimization functors
  public:
    // verbose constructor
    Linear2DFunctor(const std::vector<double>& _A,
                    const std::vector<double>& _B,
                    const std::vector<double>& _C,
                    const double _u_A, const double _u_B);
  public: 
    static std::size_t n_calls;
  
  protected:
    const std::vector<double>& A;
    const std::vector<double>& B;
    const std::vector<double>& C;
    const double u_A;
    const double u_B;
  };
  //---------------------------------------------------------------------------

  class LinearBrent2D : public Linear2DFunctor
  { // functor for brent
  public:
    // constructor
    LinearBrent2D(const std::vector<double>& _A,
                  const std::vector<double>& _B,
                  const std::vector<double>& _C,
                  const double _u_A, const double _u_B);
    
    // eval
    double operator()(double x);
  };
  //----------------------------------------------------------------------------

  class LinearNewton2D : public Linear2DFunctor
  { // functor for Newton
  public:
    // constructor
    LinearNewton2D(const std::vector<double>& _A,
                   const std::vector<double>& _B,
                   const std::vector<double>& _C,
                   const double _u_A, const double _u_B);
    
    // eval
    boost::math::tuple<double, double> operator()(double x);
  };
}

#endif // _LS_MINIMIZE_H_
