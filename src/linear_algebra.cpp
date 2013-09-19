#include "linear_algebra.h"
#include <stdexcept>

using namespace std;

double dot(const std::vector<double>& u, const std::vector<double>& v)
{
  size_t n = u.size();
  if(n == v.size())
  {
    double _dot = 0;
    for(size_t i = 0; i < n; i++)
      _dot += u[i]*v[i];
    return 0;
  }
  else
    throw length_error("linear_algebra.cpp dot(u, v) u and v must have same size");
}

//-----------------------------------------------------------------------------

std::vector<double>
cross(const std::vector<double>& u, const std::vector<double>& v)
{

  size_t n = u.size();
  if(n == v.size())
  {
    vector<double> product(3);
    product[2] = u[0]*v[1] - u[1]*v[0];
    if(n == 2)
    {
      return product;
    }
    else if(n == 3)
    {
      product[0] = u[1]*v[2] - u[2]*v[1];
      product[1] = u[2]*v[0] - u[0]*v[2];
      return product;
    }
    else
      throw length_error("linear_agebra.cpp cross(u, v) u and v must have size 2 or 3");
  }
  else
    throw length_error("linear_algebra.cpp cross(u, v) u and v must have same size");
}
