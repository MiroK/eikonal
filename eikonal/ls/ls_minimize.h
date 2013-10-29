#ifndef _LS_MINIMIZE_H_
#define _LS_MINIMIZE_H_

/*
  Local solvers based on minimizing distance + travel time functionals.
*/

#include <vector>
#include <boost/math/tools/tuple.hpp>
#include <utility>

namespace eikonal
{
  // Given triangle ABC with values of u known at A, B compute value of u at C
  // from guess u_C by monimizeing linear interpolant 
  // return to n_calls number of calls to the internal eval
  // testing interface for Brent method
  std::pair<double, double>  
  linear_brent_2d(const std::vector<double>& A,
                  const std::vector<double>& B,
                  const std::vector<double>& C,
                  const double u_A, const double u_B, const double u_C,
                  std::size_t& n_calls,
                  const std::size_t precision);
  
  // testing interface for Newton method
  std::pair<double, double>
  linear_newton_2d(const std::vector<double>& A,
                   const std::vector<double>& B,
                   const std::vector<double>& C,
                   const double u_A, const double u_B, const double u_C,
                   std::size_t& n_calls,
                   const std::size_t precision);
  
  // constructor hermite polynomial to approx current wave front
  // return the distance and grad_u_C holds the gradient at C
  std::pair<double, double>
  hermite_newton_2d(const std::vector<double>& A,
                    const std::vector<double>& B,
                    const std::vector<double>& C,
                    const double u_A, const double u_B, const double u_C,
                    const std::vector<double>& grad_u_A,
                    const std::vector<double>& grad_u_B,
                    std::vector<double>& grad_u_C,
                    std::size_t& n_calls,
                    const std::size_t precision);

  // global solver interface for Brent method
  // u_point = C, u_value is guess for the solution u_C
  // k_points = [A, B], k_values = [u[A], u[B]]
  double linear_brent_2d(const std::vector<double>& u_point,
                         const double u_value,
                         const std::vector<double>& k_points,
                         const std::vector<double>& k_values,
                         std::size_t& n_calls,
                         const std::size_t precision);
 
  // global solver interface to Newton method
  double linear_newton_2d(const std::vector<double>& u_point,
                          const double u_value,
                          const std::vector<double>& k_points,
                          const std::vector<double>& k_values,
                          std::size_t& n_calls,
                          const std::size_t precision);
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
    
    // eval df and (ddf) for we seek the minima by df = 0
    boost::math::tuple<double, double> operator()(double x);

    // return f(x)
    double eval(double x);
  };
  //---------------------------------------------------------------------------

  class HermiteNewton2D : public Linear2DFunctor
  { // functor for Newton solver of hremite...
  public:
    // constructor
    HermiteNewton2D(const std::vector<double>& _A,
                    const std::vector<double>& _B,
                    const std::vector<double>& _C,
                    const double _u_A, const double _u_B,
                    const std::vector<double>& _grad_u_A,
                    const std::vector<double>& _grad_u_B);
    
    // eval df and (ddf) for we seek the minima by df = 0
    boost::math::tuple<double, double> operator()(double x);

    // return f(x)
    double eval(double x);
    
  private:
    const std::vector<double> grad_u_A;
    const std::vector<double> grad_u_B;
  };
  //---------------------------------------------------------------------------
  
  class HermiteBrent2D : public Linear2DFunctor
  { // functor for Newton solver of hremite...
  public:
    // constructor
    HermiteBrent2D(const std::vector<double>& _A,
                    const std::vector<double>& _B,
                    const std::vector<double>& _C,
                    const double _u_A, const double _u_B,
                    const std::vector<double>& _grad_u_A,
                    const std::vector<double>& _grad_u_B);
    
    // eval
    double operator()(double x);

  private:
    const std::vector<double> grad_u_A;
    const std::vector<double> grad_u_B;
  };
  //---------------------------------------------------------------------------

  class TEST_FUNCTOR : public Linear2DFunctor
  { 
  public:
    TEST_FUNCTOR(const std::vector<double>& _A,
         const std::vector<double>& _B,
         const std::vector<double>& _C,
         const double _u_A, const double _u_B,
         const std::vector<double>& _grad_u_A,
         const std::vector<double>& _grad_u_B);
    
    boost::math::tuple<double, double> operator()(double x);

  private:
    const std::vector<double> grad_u_A;
    const std::vector<double> grad_u_B;
  };
  //---------------------------------------------------------------------------
}

#endif // _LS_MINIMIZE_H_
