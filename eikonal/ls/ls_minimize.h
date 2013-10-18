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

  // testing interface for Brent method
  double linear_brent_2d(const std::vector<double>& A,
                            const std::vector<double>& B,
                            const std::vector<double>& C,
                            const double u_A, const double u_B, const double u_C);
  
  // testing interface for Newton method
  double linear_newton_2d(const std::vector<double>& A,
                          const std::vector<double>& B,
                          const std::vector<double>& C,
                          const double u_A, const double u_B, const double u_C);
 
  // global solver interface for Brent method
  // u_point = C, u_value is guess for the solution u_C
  // k_points = [A, B], k_values = [u[A], u[B]]
  double linear_brent_2d(const std::vector<double>& u_point,
                         const double u_value,
                         const std::vector<double>& k_points,
                         const std::vector<double>& k_values);
 
  // global solver interface to Newton method
  double linear_newton_2d(const std::vector<double>& u_point,
                          const double u_value,
                          const std::vector<double>& k_points,
                          const std::vector<double>& k_values);
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
                    const double _u_A, const double _u_B) :
                        A(_A), B(_B), C(_C), u_A(_u_A), u_B(_u_B) { }

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
                  const double _u_A, const double _u_B) :
                      Linear2DFunctor(_A, _B, _C, _u_A, _u_B) { }
    
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
                   const double _u_A, const double _u_B) : 
                      Linear2DFunctor(_A, _B, _C, _u_A, _u_B) { }
    // eval
    boost::math::tuple<double, double> operator()(double x);
  };
}

#endif // _LS_MINIMIZE_H_
