from EikonalSolverLinear import EikonalSolverLinear as Linear
from EikonalSolverGeometric import EikonalSolverGeometric as Geometric
from Seeder import Segment 
from Problem import Problem
from numpy import array
from dolfin import *

if __name__ == "__main__":
  #i = 6
  #N = 2**i
  #mesh = UnitSquareMesh(N, N, "crossed")
  mesh = Mesh("../../test_results/sqr_1.msh.xml")
  #mesh = Mesh("../../test_results/mesh_perturbed_5.xml")
  
  V = FunctionSpace(mesh, "CG", 1)
  u = Function(V)
  u_exact = Function(V)
  fixed_dofs = []

  A = array([0, 0])
  B = array([1, 0])
  line_distance = Problem(Segment(A, B))
  line_distance.init(u, fixed_dofs)
  line_distance.exact_solution(u_exact)

  plot(u, interactive=True)

  solver1 = Geometric(V)
  n_sweeps = solver1.solve(u, fixed_dofs)

  l1 = assemble(abs(u - u_exact)*dx)
  l2 = assemble((u - u_exact)**2*dx)

  diff = u.vector() - u_exact.vector()
  diff.abs()
  ci = diff.max()
  
  plot(u, interactive=True)
  plot(u_exact, interactive=True)
  plot(abs(u_exact - u), interactive=True)

  print "%g %.8E %.8E %.8E %d" % (mesh.hmin(), l1, l2, ci, n_sweeps)
