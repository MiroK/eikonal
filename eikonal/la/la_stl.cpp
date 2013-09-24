#include "la_stl.h"

namespace eikonal
{
  namespace stl
  {
    double LpNorm::operator()(const double& x, const double &y)
    { 
      return x + pow(fabs(y)/_max, _p);
    }
    //-------------------------------------------------------------------------
  }
}


