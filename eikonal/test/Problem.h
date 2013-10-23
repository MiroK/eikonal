#ifndef _PROBLEM_H_
#define _PROBLEM_H_

/*
  Class for specifing problems for the Eikonal equation.
*/

#include <vector>
#include <set>
#include <string>
#include <dolfin/common/types.h>

namespace dolfin
{
  class Function;
  template<typename T> class MeshFunction;
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
    void init(std::set<dolfin::la_index>& fixed_dofs,
              dolfin::Function& u) const;

    void init(std::set<dolfin::la_index>& fixed_dofs,
              dolfin::Function& u,
              dolfin::Function& du_dx,
              dolfin::Function& du_dy) const;

    // compute the exact solution of the Eikonal equation
    void exact_solution(dolfin::Function& u) const;

    // compute the exact solution of the Eikonal equation phi, also
    // grad(u)/|grad(u) by components
    void exact_solution(dolfin::Function& u, dolfin::Function& du_dx,
                        dolfin::Function& du_dy) const;
    
    // get the CellFunction which marks cells of the mesh of u that are
    // width-cells away from the interface
    dolfin::MeshFunction<std::size_t> get_band(dolfin::Function& u,
                                       const std::size_t width) const;

    // call seeder and chose randlomly n_ref_points that are returned
    void get_ref_points(const std::size_t n_ref_points,
                        std::vector<std::vector<double> >& ref_points) const;

    // get name of the seeder
    std::string name() const;
  
  private:
    const Seeder& seeder;
  };
}

#endif // _PROBLEM_H_
