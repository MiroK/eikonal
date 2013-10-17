from dolfin import *
from numpy import dot, argmax, zeros
from math import sqrt

def accute_cell(cell, mesh):
  '''True for triangles that have no angle greater than pi/2.'''
  mesh_coordinates = mesh.coordinates()

  n = 3
  edge_sizes = zeros(3)
  vertices = cell.entities(0);
  for i in range(n):
    P = mesh_coordinates[vertices[i]]
    Q = mesh_coordinates[vertices[(i+1)%n]]
    edge_sizes[i] = sqrt(dot(P-Q, P-Q))

  i = argmax(edge_sizes)
  C = mesh_coordinates[vertices[i]]
  A = mesh_coordinates[vertices[(i+1)%n]]
  B = mesh_coordinates[vertices[(i+2)%n]]

  cosine = dot(A-C, B-C)/edge_sizes[(i+1)%n]/edge_sizes[(i+2)%n]
  return not cosine > 1

def obtuse_cells(mesh):
  '''Get indices of all cells in the mesh that are obtuse.'''
  indices = []
  for cell in cells(mesh):
    if not accute_cell(cell, mesh):
      indices.append(cell.index())
  return indices

def meshf_obtuse_cells(mesh):
  '''Return obtuse cells as mesh function.'''
  indices = obtuse_cells(mesh)

  mesh_f = CellFunction("bool", mesh)
  mesh_f.set_all(False)

  for index in indices:
    mesh_f[index] = True
  
  return mesh_f

