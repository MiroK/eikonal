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

  dofmap = V.dofmap()
  for i, x in enumerate(dofmap.tabulate_all_coordinates(mesh).reshape(V.dim(),
  2)):
    print i, x


  u = Function(V)
  u_exact = Function(V)
  fixed_dofs = []

  A = array([0, 0])
  B = array([1, 0])
  line_distance = Problem(Segment(A, B))
  line_distance.init(u, fixed_dofs)
  line_distance.exact_solution(u_exact)

  print "fixed_dofs", fixed_dofs

  plot(u, interactive=True)
  #print fixed_dofs

  solver1 = Geometric(V)
  solver1.solve(u, fixed_dofs)

  l1 = assemble(abs(u - u_exact)*dx)
  l2 = assemble((u - u_exact)**2*dx)

  plot(u, interactive=True)
  plot(u_exact, interactive=True)
  plot(abs(u_exact - u), interactive=True)

  print (u.vector() - u_exact.vector()).max()
  i_max = u.vector().array().argmax()
  print i_max
  print V.dofmap().tabulate_all_coordinates(mesh).reshape(V.dim(), 2)[i_max]
  
  print "L1 norm of error %g" % l1
  print "L2 norm of error %.16f" % l2

  #u_vector = u.vector()
  #for i, x in enumerate(u_vector):
  #  print i, "%.16f" % x


