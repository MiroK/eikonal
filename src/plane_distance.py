from dolfin import *
import numpy as np

def bottom_distance(P):
  return P[1]

def plane_distance(P, A, B, type="d"):
  t = np.dot(P - A, B - A)/np.dot(B - A, B - A)
  I = P - A - t*(B - A)
  
  if type == "sd":
    if I[0] >= 0 and I[1] >=0:
      return np.sqrt(np.dot(I, I))
    else:
      return -np.sqrt(np.dot(I, I))
  elif type == "d":
      return np.sqrt(np.dot(I, I))

def point_distance(P, point):
  return np.sqrt(np.dot(P - point, P - point))

def circle_distance(P, center, radius, type="d"):
  if type == "d":
    return abs(np.sqrt(np.dot(P-center, P-center)) - radius)
  elif type == "sd":
    return np.sqrt(np.dot(P-center, P-center)) - radius

def triangle_distance(P, A, B, C, type):
  if type == "d":
    a = plane_distance(P, A, B, "d")
    b =
    c =

#------------------------------------------------------------------------------

A = np.array([0., 1.])
B = np.array([1., 0.])
point = np.array([0.5, 0.5])
radius = 0.125

tA = np.array([1./5., 1./3.]) 
tB = np.array([4./5., 1./3.]) 
tC = np.array([0.5, 2./3.]) 

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, "CG", 2)
u = Function(V)

dofmap = V.dofmap()
for cell in cells(mesh):
  dofs = dofmap.cell_dofs(cell.index())
  dof_coordinates = dofmap.tabulate_coordinates(cell)
  for i in range(len(dofs)):
    dof = dofs[i]
    dof_x = dof_coordinates[i]

    u.vector()[dof] = triangle_distance(dof_x, A, B, C, "d")
                      #circle_distance(dof_x, point, radius, "sd")
                      #point_distance(dof_x, point)
                      #bottom_distance(dof_x)
                      #plane_distance(dof_x, A, B, "d")

plot(u, interactive=True)
