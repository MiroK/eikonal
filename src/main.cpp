#include "CG1.h"
#include "CG2.h"
#include <dolfin/generation/UnitSquareMesh.h>
#include <iostream>

#include "connectivity.h"

using namespace dolfin;

int main()
{
  UnitSquareMesh mesh(100, 100);
  CG2::FunctionSpace V(mesh);

  // build all mappings
  t_smap_la _cell_to_dof = cell_to_dof(V);
  la_vmap_t _dof_to_cell = dof_to_cell(_cell_to_dof);
  t_smap_t _cell_to_facet = cell_to_facet(V);
  t_vmap_t _facet_to_cell = facet_to_cell(_cell_to_facet);
  t_smap_t _facet_to_dof = facet_to_dof(V);
  t_vmap_t _dof_to_facet = dof_to_facet(_facet_to_dof);
  t_vmap_d _dof_to_coordinate = dof_to_coordinate(V);

  return 0;
}
