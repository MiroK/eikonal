#include "la_loop.h"
#include "la_common.h"
#include <cassert>
#include <cmath>

using namespace std;

namespace eikonal
{
  std::vector<double>
  operator+(const std::vector<double>& a, const std::vector<double>& b)
  {
    const size_t n = a.size();
    assert(n == b.size());
    vector<double> ab(a.begin(), a.end());
    
    for(size_t i = 0; i < n; i++)
      ab[i] += b[i];
    
    return ab;
  }
  //-------------------------------------------------------------------------

  std::vector<double>
  operator-(const std::vector<double>& a, const std::vector<double>& b)
  {
    const size_t n = a.size();
    assert(n == b.size());
    vector<double> ab(a.begin(), a.end());
    
    for(size_t i = 0; i < n; i++)
      ab[i] -= b[i];
    
    return ab;
  }
  //-------------------------------------------------------------------------

  std::vector<double>
  operator*(const std::vector<double>& v, const double a)
  {
    const size_t n = v.size();
    vector<double> av(v.begin(), v.end());
    
    for(size_t i = 0; i < n; i++)
      av[i] *= a;
    
    return av;
  }
  //-------------------------------------------------------------------------

  std::vector<double>
  operator*(const double a, const std::vector<double>& v)
  {
    return v*a;
  }
  //-------------------------------------------------------------------------

  std::vector<double>
  operator/(const std::vector<double>& v, const double a)
  {
    assert(a != 0.);
    return v*(1./a);
  }
  //-------------------------------------------------------------------------

  double dot(const std::vector<double>& a, const std::vector<double>& b)
  {
    const size_t n = a.size();
    assert(n == b.size());
    
    const double _max = std::max(max(abs(a)), max(abs(b)));
    double _dot = 0;
    if(_max == 0)
    {
      return _dot;
    }
    else
    {
      for(size_t i = 0; i < n; i++)
      {
        _dot += a[i]*b[i]/_max/_max;
      }
      return _dot*_max*_max;
    }
  }
  //-------------------------------------------------------------------------

  std::vector<double>
  cross(const std::vector<double>& a, const std::vector<double>& b)
  {
    assert(a.size() == b.size() and (a.size() == 2 or a.size() == 3));
    vector<double> ab(3);
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

  double norm(const std::vector<double>& v, const std::size_t p)
  {
    const size_t n = v.size();
    double _norm = 0;
    
    if(p == 0)  // maximum norm
    {
      return max(abs(v));
    }
    else
    {
      // find max element to prevent overflow
      const double max = norm(v, 0);
      if(max == 0)
      {
        return 0;
      }
      else
      {
        for(size_t i = 0; i < n; i++)
        {
          _norm += pow(fabs(v[i])/max, p);
        }
        return max*pow(_norm, 1./p);
      }
    }
  }
  //-------------------------------------------------------------------------

  double max(const std::vector<double>& v)
  {
    return v[argmax(v)];
  }
  //---------------------------------------------------------------------------

  std::size_t argmax(const std::vector<double>& v)
  {
    const size_t n = v.size();
    assert(n > 0);

    double max = v[0];
    size_t i_max = 0;
    for(size_t i = 1; i < n; i++)
    {
      double current = v[i];
      if(current > max)
      {
        i_max = i;
        max = current;
      }
    }
    return i_max;
  }
  //---------------------------------------------------------------------------

  double min(const std::vector<double>& v)
  {
    return v[argmin(v)];
  }
  //---------------------------------------------------------------------------

  std::size_t argmin(const std::vector<double>& v)
  {
    return argmax(-1*v);
  }
  //---------------------------------------------------------------------------
  
  std::vector<double> abs(const std::vector<double>& v)
  {
    std::vector<double> _v(v);
    
    for(size_t i = 0; i < v.size(); i++)
      _v[i] = abs(_v[i]);

    return _v;
  }
  //---------------------------------------------------------------------------
}
