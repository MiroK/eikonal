#ifndef _LS_MINIMIZE_H_
#define _LS_MINIMIZE_H_

/*
  Local solvers based on minimizing distance + travel time functionals.
*/

#include <vector>

namespace eikonal
{
  // Given triangle ABC with values of u known at A, B compute value of u at C
  // from guess u_C
  // uses linear interpolant to approximate position of source
  double linear_minimize_2d(const std::vector<double>& A,
                            const std::vector<double>& B,
                            const std::vector<double>& C,
                            const double u_A, const double u_B, const double u_C);
  
  // u_point = C, u_value is guess for the solution u_C
  // k_points = [A, B], k_values = [u[A], u[B]]
  double linear_minimize_2d(const std::vector<double>& u_point,
                            const double u_value,
                            const std::vector<double>& k_points,
                            const std::vector<double>& k_values);
}

namespace eikonal
{
  class Linear2DFunctor
  {
    /*
      Function to be minimized for linear interpolation on segment.
    */

    public:
      // verbos constructor
      Linear2DFunctor(const std::vector<double>& _A,
                      const std::vector<double>& _B,
                      const std::vector<double>& _C,
                      const double _u_A, const double _u_B);
      
      // eval
      double operator()(double x);

    private:
      const std::vector<double>* A;
      const std::vector<double>* B;
      const std::vector<double>* C;
      const double* u_A;
      const double* u_B;
  };
}

#endif // _LS_MINIMIZE_H_
