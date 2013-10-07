#ifndef _PROBLEM_H_
#define _PROBLEM_H_

/*
  Class for specifing problems for the Eikonal equation.
*/

#include <vector>
#include <dolfin/common/types.h>

namespace dolfin
{
  class Function;
}

namespace eikonal
{
  class Seeder;

  class Problem
  {
  public:
    // constructor
    Problem(const Seeder& _seeder);

    // seeder is used to get intersected cells, dofs in them go into
    // fixed_dofs and values of u in fixed_dofs are set to exact value
    void init(std::vector<dolfin::la_index> fixed_dofs,
              dolfin::Function& u) const;

    // compute the exact solution of the Eikonal equation
    void exact_solution(dolfin::Function& u) const;
  
  private:
    const Seeder& seeder;
  };
}

#endif // _PROBLEM_H_
