#ifndef _LOOP_H_
#define _LOOP_H_

// Loop-based functions for vector manipulation.

#include <vector>

namespace eikonal
{
  // add 2 vectors
  std::vector<double>
  operator+(const std::vector<double>& a, const std::vector<double>& b);

  // subtract 2 vectors
  std::vector<double>
  operator-(const std::vector<double>& a, const std::vector<double>& b);

  // multiply vector by scalar
  std::vector<double>
  operator*(const std::vector<double>& v, const double a);

  std::vector<double>
  operator*(const double a, const std::vector<double>& v);

  // divide vector by scalar
  std::vector<double>
  operator/(const std::vector<double>& v, const double a);

  // dot product of two vectors
  double dot(const std::vector<double>& a, const std::vector<double>& b);

  // cross product of two vectors from R^2 or R^3
  std::vector<double>
  cross(const std::vector<double>& a, const std::vector<double>& b);

  // l^p, p>-1 norm of vector, 0 is for l^\infty norm
  double norm(const std::vector<double>& v, const std::size_t p=2);

  // find the biggest element
  double max(const std::vector<double>& v);

  // find index corresponding to the biggest element
  std::size_t argmax(const std::vector<double>& v);

  // find the smallest element
  double min(const std::vector<double>& v);

  // find index corresponding to the smallest element
  std::size_t argmin(const std::vector<double>& v);
  
  // abs, x_i -> |x_i|
  std::vector<double> abs(const std::vector<double>& v);
}

#endif // _LOOP_H_
