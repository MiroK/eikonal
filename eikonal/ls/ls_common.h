#ifndef _LS_COMMON_H_
#define _LS_COMMON_H_

#include <vector>

namespace eikonal
{
  class LsData
  {
  /*
    Package of data to local solver
  */

  public:
    // constructor
    LsData(const std::vector<double>& _u_point,
           const double _u_value, 
           const std::vector<double>& _k_points, 
           const std::vector<double>& _k_values) :
                                                    u_point(_u_point),
                                                    u_value(_u_value),
                                                    k_points(_k_points),
                                                    k_values(_k_values) { }
  public:
    // coordinates of unknown dof
    const std::vector<double>& u_point;

    // guess for unknown value at unknown dof
    const double u_value;

    // coordinates of known dofs
    const std::vector<double>& k_points;

    // values at known dofs
    const std::vector<double>& k_values;
  };
}

#endif //_LS_COMMON_H_
