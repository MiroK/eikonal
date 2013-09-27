#ifndef _LS_GEOMETRIC_H_
#define _LS_GEOMETRIC_H_

/*
  Geometric solvers following "FSM for eikonal equation on tring. meshes."
*/

#include <vector>

namespace eikonal
{
  // Given triangle ABC with values of u known at A, B compute value of u at C
  // u_point = C, u_value is guess for the solution
  // k_points = [A, B], k_values = [u[A], u[B]]
  double ls_geometric_2d(const std::vector<double>& u_point,
                         const double u_value,
                         const std::vector<double>& k_points,
                         const std::vector<double>& k_values);
}

#endif // _LS_GEOMETRIC_H_
