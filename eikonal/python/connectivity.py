'''
  This module includes functionality for computing connectivity for solving
  Eikonal equation using fast sweeping method.
'''

from dolfin import cells, edges, vertices, Edge, Point

def build_cell_to_dof(V):
  '''Build a dictionary between cell and dofs that it has.'''
  mesh = V.mesh()
  dofmap = V.dofmap()
  cell_to_dof = {}
  for cell in cells(mesh):
    cell_to_dof[cell.index()] = []
    for dof in dofmap.cell_dofs(cell.index()): 
      cell_to_dof[cell.index()].append(dof)

  return cell_to_dof

def build_dof_to_cell(cell_to_dof):
  '''Invert cell_to_dof.'''
  dof_to_cell = {}
  for cell, dofs in cell_to_dof.iteritems():
    for dof in dofs:
      if dof in dof_to_cell.iterkeys():
        dof_to_cell[dof].append(cell)
      else:
        dof_to_cell[dof] = [cell]

  return dof_to_cell

def build_cell_to_edge(V):
  '''Build mapping between cell and edges that form it.'''
  cell_to_edge = {}
  mesh = V.mesh()
  mesh.init(1)
  for cell in cells(mesh):
    cell_to_edge[cell.index()] = []
    for edge in edges(cell):
      cell_to_edge[cell.index()].append(edge.index())

  return cell_to_edge

def build_edge_to_cell(cell_to_edge):
  '''Invert edge to cell.'''
  edge_to_cell = {}
  for cell, edges in cell_to_edge.iteritems():
    for edge in edges:
      if edge in edge_to_cell.iterkeys():
        edge_to_cell[edge].append(cell)
      else:
        edge_to_cell[edge] = [cell]
  return edge_to_cell

def build_dof_to_edge(edge_to_dof):
  '''Build mapping between dof and edges that have it.'''
  dof_to_edge = {}
  for edge, dofs in edge_to_dof.iteritems():
    for dof in dofs:
      if dof in dof_to_edge.iterkeys():
        dof_to_edge[dof].append(edge)
      else:
        dof_to_edge[dof] = [edge]
  return dof_to_edge

def build_edge_to_dof(V, cell_to_edge):
  '''Build mapping between edges and dofs on it.'''
  edge_to_dof = {}
  mesh = V.mesh()
  dofmap = V.dofmap()
  for cell in cells(mesh):
    cell_edges = cell_to_edge[cell.index()]
    for i in range(3): # there are 3 edges in cell
      edge = cell_edges[i]
      edge_dofs = dofmap.cell_dofs(cell.index())[dofmap.tabulate_facet_dofs(i)]
      if edge in edge_to_dof:
        for edge_dof in edge_dofs:
          edge_to_dof[edge].add(edge_dof)
      else:
        edge_to_dof[edge] = set(edge_dofs)

  return edge_to_dof

def build_dof_to_coordinate(dof_to_cell, V):
  '''Build dictionary that holds coordinates of every dof.'''
  all_coordinates = V.dofmap().tabulate_all_coordinates(V.mesh())
  gdim = V.mesh().geometry().dim()
  dof_to_coordinate = {}
  for i, dof in enumerate(dof_to_cell.keys()):
    dof_to_coordinate[dof] = all_coordinates[i*gdim : (i+1)*gdim]

  return dof_to_coordinate

def build_dof_to_dof(cell_to_dof, dof_to_cell):
  '''Build connectivity map between degress of freedom.'''
  dof_to_dof = {}
  for dof in dof_to_cell.iterkeys():
    dof_to_dof[dof] = set([dof])
    for cell in dof_to_cell[dof]:
      for _dof in cell_to_dof[cell]:
        dof_to_dof[dof].add(_dof)
    dof_to_dof[dof].remove(dof)

  return dof_to_dof


if __name__ == '__main__':
  from dolfin import *
  mesh = UnitSquareMesh(2, 2)
  V = FunctionSpace(mesh, "CG", 1)
  cell_dof = build_cell_to_dof(V)
  dof_cell = build_dof_to_cell(cell_dof)
  dof_dof = build_dof_to_dof(cell_dof, dof_cell)

  plot(mesh, interactive=True)
  print "dof->coordinate", build_dof_to_coordinate(dof_cell, V)
  print
  print "cell->dof: ", cell_dof
  print
  print "dof->cell: ", dof_cell
  print
  print "dof->dof: ", dof_dof

