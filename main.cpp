#include "eikonal.h"
#include <cstdlib>
#include <dolfin.h>
#include <iostream>
#include "test.h"

using namespace eikonal;
using namespace dolfin;

class MyExpression : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 1;
  }
};

int main(int argc, char* argv[])
{
  double _c1[2] = {-1., 0.}; std::vector<double> c1(_c1, _c1+2);
  double _c2[2] = {sqrt(1.5), 0}; std::vector<double> c2(_c2, _c2+2);
  TwoCircles two_circles("twocircle", c1, 0.5, c2, 0.5);
  Problem problem(two_circles);

  RectangleMesh mesh(-2, -2, 2, 2, 40, 40);
  test::CoefficientSpace_u V(mesh);
  Function u(V);
  MeshFunction<std::size_t> band = problem.get_band(u, 3);
  plot(band);
  interactive(true);


  /*
  UnitSquareMesh mesh(10, 10);
  std::size_t dim = mesh.topology().dim();
  MeshFunction<std::size_t> foo(mesh, dim);
  foo.set_all(0);
  foo[110] = 1.;

  mesh.init(dim, dim);
  for(std::size_t level = 0; level < 4; level++)
  {
    plot(foo);
    interactive(true);
    std::set<std::size_t> marked_cells;
    for(CellIterator cell(mesh); !cell.end(); ++cell)
    {
      if(foo[*cell] == 1)
      {
        std::size_t num_next = cell->num_entities(dim);
        const std::size_t* next = cell->entities(dim);
        marked_cells.insert(next, next+num_next);
      }
    }
    std::set<std::size_t>::const_iterator marked_cell = marked_cells.begin();
    for(; marked_cell != marked_cells.end(); marked_cell++)
    {
      foo[*marked_cell] = 1;
    }
  }

  SubMesh sub_mesh(mesh, foo, 1);
  plot(sub_mesh);
  interactive(true);
  
  Constant one(1.);
  test::Form_area area(sub_mesh, one);
  double a = assemble(area);
  std::cout << a << std::endl;

  test::CoefficientSpace_u V(mesh);
  Function u(V);
  u.interpolate(MyExpression());
  plot(u);
  interactive(true);

  test::Form_norm norm(sub_mesh, u);
  double n = assemble(norm);
  std::cout << n << std::endl;
  */
  //all_hermite_tests(2);  
  //return 0;
}

