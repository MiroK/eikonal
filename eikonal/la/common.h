#ifndef _COMMON_H_
#define _COMMON_H_

// Common functionality for vectors.

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

namespace eikonal
{
  // equality tolerance for comparison
  const double LA_EPS = 1E-12;
  
  // comparison of two entries
  template<typename T>
  bool compare(const T& a, const T& b, const double EPS=LA_EPS);

  // templeted absolute value
  template<typename T> T abs(const T& a);

  // print vector
  template<typename T>
  void print(const std::vector<T>& v);

  // comparison, uses compare function with the tolarance
  template<typename T>
  bool operator==(const std::vector<T>& a, const std::vector<T>& b);
}

// definitions

namespace eikonal
{
  template<typename T>
  bool compare(const T& a, const T& b, const double EPS)
  {
    return abs(a-b) < EPS;
  }
  //---------------------------------------------------------------------------

  template<typename T>
  void print(const std::vector<T>& v)
  {
    std::cout << "[";
    for(size_t i = 0; i < v.size()-1; i++)
      std::cout << v[i] << ", ";

    std::cout << v[v.size()-1] <<"]\n";
  }
  //---------------------------------------------------------------------------

  template<typename T> T abs(const T& a)
  {
    return a > 0 ? a : -a; 
  }
  //---------------------------------------------------------------------------

  template<typename T>
  bool operator==(const std::vector<T>& a, const std::vector<T>& b)
  {
    assert(a.size() == b.size()); // no stl here
    for(size_t i = 0; i < a.size(); i++)
    {
      if(not compare(a[i], b[i]))
        return false;
    }
    return true;
  }
  //---------------------------------------------------------------------------
}


#endif // _COMMON_H_
