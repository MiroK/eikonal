#ifndef _LINEAR_ALGEBRA_H_
#define _LINEAR_ALGEBRA_H_

/*
  Some functions from linear algebra.
*/

#include <vector>

// dot product of two vectors of arbitrary size
double dot(const std::vector<double>& u, const std::vector<double>& v);

// cross product of vectors of size 2, 3
std::vector<double>
cross(const std::vector<double>& u, const std::vector<double>& v);

#endif //_LINEAR_ALGEBRA_H_
