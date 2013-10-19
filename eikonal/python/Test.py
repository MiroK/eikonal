from EikonalSolverLinear import EikonalSolverLinear as Linear
from EikonalSolverQuadratic import EikonalSolverQuadratic as Quadratic
from EikonalSolverGeometric import EikonalSolverGeometric as Geometric
from EikonalSolverSampleByLine import EikonalSolverSampleByLine as SampleLine
from Seeder import Segment, TwoCircle, MyPoint
from Problem import Problem
from accute_statistics import obtuse_cells
from numpy import array
from math import sqrt
from time import clock
from dolfin import *

# 0 -- Rectangle right, 1 -- Rectangle crossed, 2 -- Gmsh
def two_circle(i, mesh_type, xtol, plot_on=False):
  names = ["rect_right", "rect_crossed", "rect_perturbed", "rect_gmsh"]

  print "two_circle", names[mesh_type]

  if mesh_type == 0:
    N = 2**i
    mesh = RectangleMesh(-2, -2, 2, 2, N, N)
  elif mesh_type == 1:
    N = 2**i
    mesh = RectangleMesh(-2, -2, 2, 2, N, N, 'crossed')
  elif mesh_type == 2:
    mesh = Mesh("../../test_results/perturbed_rect_%d.xml" % i)
  elif mesh_type == 3:
    mesh = Mesh("../../test_results/rectangle_%d.msh.xml" % i)

  n_cells = mesh.num_cells()
  n_obtuse_cells = len(obtuse_cells(mesh))

  V = FunctionSpace(mesh, "CG", 1)
  u = Function(V)
  u_exact = Function(V)
  fixed_dofs = []

  c1 = array([-1., 0.])
  c2 = array([sqrt(1.5), 0.])
  r = 0.5
  two_circle_distance = Problem(TwoCircle(c1, r, c2, r))
  two_circle_distance.init(u, fixed_dofs)
  two_circle_distance.exact_solution(u_exact)

  if plot_on:
    plot(u_exact, interactive=True)
    plot(u, interactive=True)

  #solver1 = Geometric(V)
  solver1 = Linear(V)
  start = clock()
  n_sweeps = solver1.solve(u, fixed_dofs, xtol=xtol)
  stop = clock() - start

  l1 = assemble(abs(u - u_exact)*dx)
  l2 = assemble((u - u_exact)**2*dx)

  diff = u.vector() - u_exact.vector()
  diff.abs()
  ci = diff.max()
  
  if plot_on:
    plot(u, interactive=True)
    plot(u_exact, interactive=True)
    plot(abs(u_exact - u), interactive=True)

  if hasattr(solver1, "max_calls") and hasattr(solver1, "min_calls"):
    print "%g %.8E %.8E %.8E %d/%d %d %g %d %d" %\
    (mesh.hmin(), l1, l2, ci, n_obtuse_cells, n_cells, n_sweeps, stop,
     solver1.min_calls, solver1.max_calls)
  else:
    print "%g %.8E %.8E %.8E %d/%d %d %g" %\
    (mesh.hmin(), l1, l2, ci, n_obtuse_cells, n_cells, n_sweeps, stop)

#-------------------------------------------------------------------------------

# 0 -- Rectangle right, 1 -- Rectangle crossed, 2 -- Gmsh
def line(i, mesh_type, xtol, plot_on=False):
  names = ["rect_right", "rect_crossed", "rect_perturbed", "rect_gmsh"]

  if mesh_type == 0:
    N = 2**i
    mesh = UnitSquareMesh(N, N)
    #mesh = RectangleMesh(-2, -2, 2, 2, N, N)
  #elif mesh_type == 1:
  #  N = 2**i
  #  mesh = RectangleMesh(-2, -2, 2, 2, N, N, 'crossed')
  #elif mesh_type == 2:
  #  mesh = Mesh("../../test_results/perturbed_rect_%d.xml" % i)
  #elif mesh_type == 3:
  #  mesh = Mesh("../../test_results/rectangle_%d.msh.xml" % i)

  n_cells = mesh.num_cells()
  n_obtuse_cells = len(obtuse_cells(mesh))
  
  V = FunctionSpace(mesh, "CG", 2)
  
  u = Function(V)
  u_exact = Function(V)
  fixed_dofs = []

  #
  A = array([1., 0.])
  point_distance = Problem(MyPoint(A))
  point_distance.init(u, fixed_dofs)
  point_distance.exact_solution(u_exact)

  # distance from line
  #A = array([0, 0.1])
  #B = array([1, 0.1])
  #line_distance = Problem(Segment(A, B))
  #line_distance.init(u, fixed_dofs)
  #line_distance.exact_solution(u_exact)

  if plot_on:
    plot(u_exact, interactive=True)
    plot(u, interactive=True)

  #print fixed_dofs

  #solver1 = Geometric(V)
  #solver1 = Linear(V)
  solver1 = Quadratic(V)
  #solver1 = SampleLine(V)

  start = clock()
  n_sweeps = solver1.solve(u, fixed_dofs, xtol=xtol)
  stop = clock() - start

  l1 = assemble(abs(u - u_exact)*dx)
  l2 = assemble((u - u_exact)**2*dx)

  diff = u.vector() - u_exact.vector()
  #print u.vector().array()
  diff.abs()
  ci = diff.max()
  
  if plot_on:
    plot(u, interactive=True)
    plot(u_exact, interactive=True)
    plot(abs(u_exact - u), interactive=True)

  if hasattr(solver1, "max_calls") and hasattr(solver1, "max_calls"):
    print "%g %.8E %.8E %.8E %d/%d %d %g %d %d" %\
    (mesh.hmin(), l1, l2, ci, n_obtuse_cells, n_cells, n_sweeps, stop,
     solver1.min_calls, solver1.max_calls)
  else:
    print "%g %.8E %.8E %.8E %d/%d %d %g" %\
    (mesh.hmin(), l1, l2, ci, n_obtuse_cells, n_cells, n_sweeps, stop)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
  for i in range(1, 6):
    line(i, 0, 1E-12, False)
#  for i in range(3, 8):
#    two_circle(i, 1, False)
#  print
#
#  for i in range(3, 8):
#    two_circle(i, 2, False)
#  print
#
#  for i in range(1, 6):
#    two_circle(i, 3, False)
#  print

  

