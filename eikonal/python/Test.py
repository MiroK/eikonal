from EikonalSolverLinear import EikonalSolverLinear as Linear
from EikonalSolverGeometric import EikonalSolverGeometric as Geometric
from Seeder import Segment 
from Problem import Problem
from numpy import array
from dolfin import *

if __name__ == "__main__":
  mesh = UnitSquareMesh(2, 2)
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
  #print fixed_dofs

  solver1 = Geometric(V)
  solver1.solve(u, fixed_dofs)

  l1 = assemble(abs(u - u_exact)*dx)
  l2 = assemble((u - u_exact)**2*dx)

  plot(u, interactive=True)
  plot(u_exact, interactive=True)

  print "L1 norm of error", l1
  print "L2 norm of error", l2
