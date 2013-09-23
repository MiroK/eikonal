#ifndef _STL_H_
#define _STL_H_

// STL-based functions for vector manipulation.

#include <vector>
#include <cassert>
#include <algorithm>   // transform
#include <functional>  // plus, minus, multiplies, bind1st
#include <numeric>     // inner_product
#include "common.h" // abs

namespace eikonal
{
  namespace stl
  {
    // add 2 vectors
    template<typename T> std::vector<T>
    operator+(const std::vector<T>& a, const std::vector<T>& b);

    // subtract 2 vectors
    template<typename T> std::vector<T>
    operator-(const std::vector<T>& a, const std::vector<T>& b);

    // multiply vector by scalar
    template<typename T> std::vector<T>
    operator*(const std::vector<T>& v, const double a);

    template<typename T> std::vector<T>
    operator*(const double a, const std::vector<T>& v);

    // divide vector by scalar
    template<typename T> std::vector<T>
    operator/(const std::vector<T>& v, const double a);

    // dot product of two vectors
    template<typename T> double
    dot(const std::vector<T>& a, const std::vector<T>& b);

    // cross product of two vectors from R^2 or R^3
    template<typename T> std::vector<T>
    cross(const std::vector<T>& a, const std::vector<T>& b);

    // l^p, p>-1 norm of vector, 0 is for l^\infty norm
    template<typename T> double
    norm(const std::vector<T>& v, const std::size_t p=2);

    // norm computation hack
    class LpNorm
    {
    public:
      // constuctor, maximum element of vector, norm exponent
      LpNorm(const double max, const std::size_t p) : _max(max), _p(p) { }

      // accumulator
      double operator()(const double& x, const double &y);

    private:
      const double _max;
      const std::size_t _p;
    };

    // get biggest element of vector
    template<typename T> T max(const std::vector<T>& v);
    
    // get index of the biggest element
    template<typename T> std::size_t argmax(const std::vector<T>& v);

    // get smallest element of vector
    template<typename T> T min(const std::vector<T>& v);

    // get index of the smallest element
    template<typename T> std::size_t argmin(const std::vector<T>& v);

    // apply absolute value to all elements of the vector
    template<typename T> std::vector<T> abs(const std::vector<T>& v);
  }
}

// definitions

namespace eikonal
{
  namespace stl
  {
    template<typename T> std::vector<T>
    operator+(const std::vector<T>& a, const std::vector<T>& b)
    {
      assert(a.size() == b.size());
      std::vector<T> ab(a.begin(), a.end());
      std::transform(b.begin(), b.end(), ab.begin(), ab.begin(), std::plus<T>());
      return ab;
    }
    //-------------------------------------------------------------------------

    template<typename T> std::vector<T>
    operator-(const std::vector<T>& a, const std::vector<T>& b)
    {
      assert(a.size() == b.size());
      std::vector<T> ab(a.begin(), a.end());
      std::transform(ab.begin(), ab.end(), b.begin(), ab.begin(), std::minus<T>());
      return ab;
    }
    //-------------------------------------------------------------------------

    template<typename T> std::vector<T>
    operator*(const std::vector<T>& v, const double a)
    {
      std::vector<T> av(v.begin(), v.end());
      std::transform(av.begin(), av.end(), av.begin(),
                std::bind1st(std::multiplies<T>(), a));
      return av;
    }
    //-------------------------------------------------------------------------

    template<typename T> std::vector<T>
    operator*(const double a, const std::vector<T>& v)
    {
      return v*a;
    }
    //-------------------------------------------------------------------------

    template<typename T> std::vector<T>
    operator/(const std::vector<T>& v, const double a)
    {
      assert(a != 0);
      return v*(1./a);
    }
    //-------------------------------------------------------------------------

    template<typename T> double
    dot(const std::vector<T>& a, const std::vector<T>& b)
    {
      assert(a.size() == b.size());
      
      T max_ab = std::max(max(abs(a)), max(abs(b)));
      if(max_ab == 0)
      {
        return 0;
      }
      else
      { // it is possible to save memory by custom accumulator
        std::vector<T> _a = a/max_ab;
        std::vector<T> _b = b/max_ab;
      
        return std::inner_product(_a.begin(),
                                  _a.end(), _b.begin(), 0.0)*max_ab*max_ab;
      }
    }
    //-------------------------------------------------------------------------

    template<typename T> std::vector<T>
    cross(const std::vector<T>& a, const std::vector<T>& b)
    {
      assert(a.size() == b.size() and (a.size() == 2 or a.size() == 3));
      std::vector<T> ab(3);
      ab[2] = a[0]*b[1] - a[1]*b[0];
      if(a.size() == 2)
        return ab;
      else
      {
        ab[0] = a[1]*b[2] - a[2]*b[1];
        ab[1] = a[2]*b[0] - a[0]*b[2];
        return ab;
      }
    }
    //-------------------------------------------------------------------------

    template<typename T> double
    norm(const std::vector<T>& v, const std::size_t p)
    {
      if(p == 0)
      {
        return std::max(*std::max_element(v.begin(), v.end()),
                        eikonal::abs(*std::min_element(v.begin(), v.end())));
      }
      else
      {
        const double max = norm(v, 0);
        if(max == 0)
        {
          return 0;
        }
        else
        {
          LpNorm lp_norm(max, p);
          double sum = std::accumulate(v.begin(), v.end(), 0.0, lp_norm);
          return max*pow(sum, 1./p);
        }
      }
    }
    //-------------------------------------------------------------------------

    template<typename T> T max(const std::vector<T>& v)
    {
      return *std::max_element(v.begin(), v.end());
    }
    //-------------------------------------------------------------------------
    
    template<typename T> std::size_t argmax(const std::vector<T>& v)
    {
      return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
    }
    //-------------------------------------------------------------------------

    template<typename T> T min(const std::vector<T>& v)
    {
      return *std::min_element(v.begin(), v.end());
    }
    //-------------------------------------------------------------------------

    template<typename T> std::size_t argmin(const std::vector<T>& v)
    {
      return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
    }
    //-------------------------------------------------------------------------

    template<typename T> std::vector<T> abs(const std::vector<T>& v)
    {
      std::vector<T> _v(v);
      std::transform(_v.begin(), _v.end(), _v.begin(), eikonal::abs<T>);
      return _v;
    }
    //-------------------------------------------------------------------------
  }
}

#endif // _STL_H_
