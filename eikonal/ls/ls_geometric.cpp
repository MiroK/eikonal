#include "ls_geometric.h"
#include "la/la_common.h"
#include "la/la_loop.h"
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;

namespace eikonal
{
  double linear_geometric_2d(const std::vector<double>& A,
                             const std::vector<double>& B,
                             const std::vector<double>& C,
                             const double u_A, const double u_B, const double u_C)
  {
    const double a = norm(C - B);
    const double b = norm(C - A);
    const double c = norm(A - B);

    if(abs(u_B - u_A) <= c)
    {
      const double theta = asin(abs(u_B - u_A)/c);
      const double alpha = acos(dot(C - B, A - B)/a/c);
      const double beta = acos(dot(C - A, B - A)/b/c);
      const double pip = M_PI/2.;

      if((std::max(0., alpha - pip) <= theta and theta <= (pip - beta)) or
         ((alpha - pip) <= theta and theta <= std::min(0., pip - beta)))
      {
        const double h = a*sin(alpha - theta);
        const double H = b*sin(beta + theta);
        
        vector<double> values(2);
        values[0] = u_C;
        values[1] = 0.5*(h + u_B) + 0.5*(H + u_A);
        return *min_element(values.begin(), values.end());
      }
      else
      {
        vector<double> values(3);
        values[0] = u_C;
        values[1] = u_A + b;
        values[2] = u_B + a;
        return *min_element(values.begin(), values.end());
      }
    }
    else
    {
      vector<double> values(3);
      values[0] = u_C;
      values[1] = u_A + b;
      values[2] = u_B + a;
      return *min_element(values.begin(), values.end());
    }
  }
  //---------------------------------------------------------------------------
  
  double linear_geometric_2d(const std::vector<double>& u_point,
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
    return linear_geometric_2d(A, B, C, u_A, u_B, u_C);
  }
  //---------------------------------------------------------------------------
}
