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
    const double a = norm(C - B, 2);
    const double b = norm(C - A, 2);
    const double c = norm(A - B, 2);

    if((u_B - u_A) <= c)
    {
      double theta = asin((u_B - u_A)/c);
      const double alpha = acos(dot(C - B, A - B)/a/c);
      const double beta = acos(dot(C - A, B - A)/b/c);
      const double pip = M_PI/2.;

      std::cout << "\tu_A " << u_A << std::endl; 
      std::cout << "\tu_B " << u_B << std::endl;
      std::cout << "\tA "; print(A);
      std::cout << "\tB "; print(B);
      std::cout << "\talpha " << alpha << std::endl;
      std::cout << "\tbeta " << beta << std::endl;
      std::cout << "\ttheta " << theta << std::endl;
      std::cout << "\ta " << a << std::endl;
      std::cout << "\tb " << b << std::endl;
      std::cout << "\tc " << c << std::endl;
      
      if((std::max(0., alpha - pip) <= theta and theta <= (pip - beta)))
      {
        const double h = a*sin(alpha - theta);

        vector<double> values(2);
        values[0] = u_C;
        values[1] = h + u_B;

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
    const std::size_t i = k_values[0] <= k_values[1] ? 0 : 1;
    const std::size_t ip1 = (i+1)%2;
    const double u_A = k_values[i], u_B = k_values[ip1], u_C = u_value;
    
    const std::vector<double> A(&k_points[i*dim], &k_points[i*dim] + dim),
                              B(&k_points[ip1*dim], &k_points[ip1*dim] + dim),
                              C(u_point);
    
    // use the verbose interface
    return linear_geometric_2d(A, B, C, u_A, u_B, u_C);
  }
  //---------------------------------------------------------------------------
}
