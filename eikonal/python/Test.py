from EikonalSolverLinear import EikonalSolverLinear as Linear
from EikonalSolverQuadratic import EikonalSolverQuadratic as Quadratic
from EikonalSolverGeometric import EikonalSolverGeometric as Geometric
from Quad import Quad
from QuadraticRatic import QuadraticRatic
from Q4 import Q4

from EikonalSolverSampleByLine import EikonalSolverSampleByLine as SampleLine
from Seeder import Segment, TwoCircle, MyPoint
from Problem import Problem
from accute_statistics import obtuse_cells
from numpy import array
from math import sqrt
from time import clock
from dolfin import *

def test(i, solver, choice, problem):
  xtol = 1E-12
  N = 2**i
  mesh = UnitSquareMesh(N, N, 'crossed')
  
  n_cells = mesh.num_cells()
  n_obtuse_cells = len(obtuse_cells(mesh))

  if problem == "point":
    A = array([1.0, 0.0])
    point_distance = Problem(MyPoint(A))
  elif problem == "line":
    A = array([0., 0.])
    B = array([1., 0.])
    point_distance = Problem(Segment(A, B))

  if solver == "Q1":
    # set up the linear problem
    V = FunctionSpace(mesh, "CG", 1)
    u = Function(V)
    u_exact = Function(V)
    fixed_dofs = []
    point_distance.init(u, fixed_dofs)
    point_distance.exact_solution(u_exact)

    solver1 = Geometric(V)
    start = clock()
    n_sweeps = solver1.solve(u, fixed_dofs, xtol=xtol)
    stop = clock() - start
    l1 = assemble(abs(u - u_exact)*dx)
    l2 = assemble((u - u_exact)**2*dx)

    # set up the quadratic problem
    W = FunctionSpace(mesh, "CG", 2)
    v = interpolate(u, W)
    v_exact = Function(W)
    fixed_dofs = []
    point_distance.init(v, fixed_dofs)
    point_distance.exact_solution(v_exact)

    solver2 = Quadratic(W, choice)
    start = clock()
    n_sweeps = solver2.solve(v, fixed_dofs, xtol=xtol, order_vector=v.vector())
    stop = clock() - start

  if solver == "Q2":
    # set up the linear problem
    V = FunctionSpace(mesh, "CG", 1)
    u = Function(V)
    u_exact = Function(V)
    fixed_dofs = []
    point_distance.init(u, fixed_dofs)
    point_distance.exact_solution(u_exact)

    solver1 = Geometric(V)
    n_sweeps = solver1.solve(u, fixed_dofs, xtol=xtol)

    # set up the quadratic problem
    W = FunctionSpace(mesh, "CG", 2)
    v = interpolate(u, W)
    xxx = Function(W)
    xxx.vector()[:] = v.vector().array()

    v_exact = Function(W)
    Fixed_dofs = []
    point_distance.init(v, Fixed_dofs)
    point_distance.exact_solution(v_exact)

    #print fixed_dofs
    #for i in range(W.dim()):
    #  print i, xxx.vector()[i], v.vector()[i],
    #  if i in fixed_dofs:
    #    print "<--"
    #  else:
    #    print



    solver2 = Quad(W, choice)
    start = clock()
    n_sweeps = solver2.solve(v, Fixed_dofs, xtol=xtol)
    stop = clock() - start
  
  if solver == "Q3":
    # set up the quadratic problem
    W = FunctionSpace(mesh, "CG", 2)
    v = Function(W) 
    v_exact = Function(W)
    fixed_dofs = []
    point_distance.init(v, fixed_dofs)
    point_distance.exact_solution(v_exact)

    solver2 = QuadraticRatic(W, choice)
    start = clock()
    n_sweeps = solver2.solve(v, fixed_dofs, xtol=xtol)
    stop = clock() - start
  
  if solver == "Q4":
    # set up the linear problem
    V = FunctionSpace(mesh, "CG", 1)
    u = Function(V)
    u_exact = Function(V)
    fixed_dofs = []
    point_distance.init(u, fixed_dofs)
    point_distance.exact_solution(u_exact)

    solver1 = Geometric(V)
    start = clock()
    n_sweeps = solver1.solve(u, fixed_dofs, xtol=xtol)
    stop = clock() - start
    l1 = assemble(abs(u - u_exact)*dx)
    l2 = assemble((u - u_exact)**2*dx)

    # set up the quadratic problem
    W = FunctionSpace(mesh, "CG", 2)
    v = Function(W) 
    v_exact = Function(W)
    fixed_dofs = []
    point_distance.init(v, fixed_dofs)
    point_distance.exact_solution(v_exact)
    order = interpolate(u, W)

    solver2 = Q4(W, choice)
    start = clock()
    n_sweeps = solver2.solve(v, fixed_dofs, xtol=xtol, order_u=order)
    stop = clock() - start
 
  l1 = assemble(abs(v - v_exact)*dx)
  l2 = assemble((v - v_exact)**2*dx)
  diff = v.vector() - v_exact.vector()
  diff.abs()
  ci = diff.max()
  
  plot(v, interactive=True, title='numeric')
  plot(v_exact, interactive=True, title='exact')
  
  v_exact.vector()[:] -= v.vector().array()
  v_exact.vector().abs()
  plot(v_exact, interactive=True, title='error')
 
  if choice == 0:
    local = "sample"
  else:
    local = "quad"

  print "%g %.8E %.8E %.8E %d/%d %d %g" %\
    (mesh.hmin(), l1, l2, ci, n_obtuse_cells, n_cells, n_sweeps, stop)

  out_name = "./test_results/"+solver+"_"+local+"_"+problem+".dat"
  with open(out_name, "a") as out:
    out.write("%g %.8E %.8E %.8E %d %d %d %g\n" %\
    (mesh.hmin(), l1, l2, ci, n_obtuse_cells, n_cells, n_sweeps, stop)
    )
#-------------------------------------------------------------------------------

if __name__ == "__main__":
  for solver in ['Q2']:
    print solver
    for choice in [1]:
      for i in range(1, 3):
        test(i, solver=solver, choice=choice, problem="point")
      print
  

