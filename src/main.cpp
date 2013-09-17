#include "CG1.h"
#include <dolfin/generation/UnitSquareMesh.h>

#include "connectivity.h"

using namespace dolfin;

int main()
{
  UnitSquareMesh mesh(2, 2);
  CG1::FunctionSpace V(mesh);

  // build cell_to_dof map and see if it is okay
  tMapLa _cell_to_dof = cell_to_dof(V);
  info("%s", print<tMapLa, la_index>(_cell_to_dof, "Cell->Dof").c_str());

  laMapT test = invert_map<tMapLa, laMapT, la_index>(_cell_to_dof);

  // build dof_to_cell map and see if it is okay
  laMapT _dof_to_cell = dof_to_cell(_cell_to_dof);
  info("%s", print<laMapT, std::size_t>(_dof_to_cell, "Dof->Cell").c_str());


  info("%s", print<laMapT, std::size_t>(test, "Dof->Cell1").c_str());

  // build cell_to_facet map and see if it is okay
  tMapT _cell_to_facet = cell_to_facet(V);
  info("%s", print<tMapT, std::size_t>(_cell_to_facet, "Cell->Edge").c_str());

  // build facet_to_cell and see if it is okay
  tMapT _facet_to_cell = facet_to_cell(_cell_to_facet);
  info("%s", print<tMapT, std::size_t>(_facet_to_cell, "Edge->Cell").c_str());
 
  // build facet_to_dof and see if it is okay
  tMapT _facet_to_dof = facet_to_dof(V);
  info("%s", print<tMapT, std::size_t>(_facet_to_dof, "Edge->Dof").c_str());

  // build dof_to_facet and see if it is okay
  tMapT _dof_to_facet = dof_to_facet(_facet_to_dof);
  info("%s", print<tMapT, std::size_t>(_dof_to_facet, "Dof->Edge").c_str());
  
  return 0;
}
