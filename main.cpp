#include <dolfin/generation/UnitSquareMesh.h>
#include "eikonal.h"

using namespace dolfin;
using namespace eikonal;

int main()
{
  UnitSquareMesh mesh(2, 2);
  CG1::FunctionSpace V(mesh);

  Solver solver(V);

  bool _ds[9] = {true, true, false, false, true, true, false, false, true};
  std::vector<bool> dof_status(_ds, _ds + 9);
  //print(solver.get_cell_set_dofs(5, 8, dof_status));

  std::vector<std::vector<double> > ref_points;
  std::vector<double> p0; p0.push_back(1.); p0.push_back(0.); 
  std::vector<double> p1; p1.push_back(0.5); p1.push_back(0.5); 
  std::vector<double> p2; p2.push_back(0.); p2.push_back(0.); 
  ref_points.push_back(p0);
  ref_points.push_back(p1);
  ref_points.push_back(p2);
  
  la_index _dofs[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<la_index> dofs(_dofs, _dofs + 9);

  t_vmap_d map = dof_to_coordinate(V);
  std::cout << print_map<t_vmap_d, std::vector<double> >(map, "x"); 

  LpDistanceSorter sorter(map);
  sorter.sort(dofs, ref_points, 2);

  for(std::size_t k = 0; k < ref_points.size(); k++)
  {
    std::cout << "Dofs sorted by distance from ";
    print(ref_points[k]);
    for(std::size_t reverse = 0; reverse < 2; reverse++)
    {
      std::cout << " in " << reverse << " order: ";
      for(MyIterator<la_index> dof = sorter.get(k, reverse); !dof.end(); ++dof)
      {
        std::cout << *dof << " ";
      }
      std::cout << std::endl;
    }
  }
  
  return 0;
}
