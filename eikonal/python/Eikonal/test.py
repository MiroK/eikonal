'''
  Convergence rate test for eikonal solvers that use fast sweeping method.
'''

from EikonalSolverLinear import EikonalSolverLinear as Linear
from EikonalSolverQuadratic import EikonalSolverQuadratic as Quadratic
import numpy as np
import sys
from problem import *

def test_rate_linear():
  '''
    Convergence rate of the linear Eikonal solver. Fenics and Gmsh meshes are
    used.
  '''
  for k in range(2):
    if k == 0:
      meshes = [Mesh("meshes/rect_%d.msh.xml" % i) for i in [0, 1, 2, 3]]
      f = open("linear_rates_gmsh.txt", "w")
    else:
      meshes = [RectangleMesh(-2, -2, 2, 2, N, N) for N in\
                                          [8, 16, 32, 64, 128, 256]]
      f = open("linear_rates_fenics.txt", "w")

    hs = [] # mesh size
    iters = [] # number of iterations
    errors0 = [] # linfty errors
    errors1 = [] # l1 errors
    errors2 = [] # errornorm error

    for mesh in meshes:
      h, errors, num_iters = two_circle_problem(mesh, 1, Linear)
      hs.append(h)
      errors0.append(errors[0])
      errors1.append(errors[1])
      errors2.append(errors[2])
      iters.append(num_iters)

    for i in range(1, len(meshes)):
      r0 = float(numpy.log(errors0[i]/errors0[i-1])/numpy.log(hs[i]/hs[i-1]))
      r1 = float(numpy.log(errors1[i]/errors1[i-1])/numpy.log(hs[i]/hs[i-1]))
      r2 = float(numpy.log(errors2[i]/errors2[i-1])/numpy.log(hs[i]/hs[i-1]))
      num_iters = iters[i]
      f.write("l_infty=%g\tl1=%g\tl2=%g\t\t(%d)\n" % (r0, r1, r2, num_iters))
    f.close()

def test_rate_quadratic():
  '''
    Convergence rate of the quadratic Eikonal solver. Fenics and Gmsh meshes are
    used.
  '''
  for k in range(2):
    if k == 0:
      meshes = [Mesh("meshes/rect_%d.msh.xml" % i) for i in [0, 1, 2, 3]]
      f = open("quadratic_rates_gmsh.txt", "w")
    else:
      meshes = [RectangleMesh(-2, -2, 2, 2, N, N) for N in\
                                          [8, 16, 32, 64, 128, 256]]
      f = open("quadratic_rates_fenics.txt", "w")

    hs = [] # mesh size
    iters = [] # number of iterations
    errors0 = [] # linfty errors
    errors1 = [] # l1 errors
    errors2 = [] # errornorm error

    for mesh in meshes:
      h, errors, num_iters = two_circle_problem(mesh, 2, Quadratic)
      hs.append(h)
      errors0.append(errors[0])
      errors1.append(errors[1])
      errors2.append(errors[2])
      iters.append(num_iters)

    for i in range(1, len(meshes)):
      r0 = float(numpy.log(errors0[i]/errors0[i-1])/numpy.log(hs[i]/hs[i-1]))
      r1 = float(numpy.log(errors1[i]/errors1[i-1])/numpy.log(hs[i]/hs[i-1]))
      r2 = float(numpy.log(errors2[i]/errors2[i-1])/numpy.log(hs[i]/hs[i-1]))
      num_iters = iters[i]
      f.write("l_infty=%g\tl1=%g\tl2=%g\t\t(%d)\n" % (r0, r1, r2, num_iters))
    f.close()

#-----------------------------------------------------------------------------

if __name__ == "__main__":
  if sys.argv[1] == "linear":
    test_rate_linear()
  elif sys.argv[1] == "quadratic":
    test_rate_quadratic()

#mesh = UnitSquareMesh(3, 3)
#V = FunctionSpace(mesh, "CG", 2)
#solver = Quadratic(V)

#print "cell to dof", solver.cell_to_dof
#print "dof to cell", solver.dof_to_cell
#print "edge to cell", solver.edge_to_cell
#print "cell to edge", solver.cell_to_edge
#print "edge to dof:", solver.edge_to_dof
#print "dof to edge:", solver.dof_to_edge
#print "dot to coordinate", solver.dof_to_coordinate
#u = Function(V)
#u.vector()[:] = 100
#u.vector()[0] = 0
#u.vector()[2] = 0
#u.vector()[5] = 0


#fixed_dofs = [0, 2, 5]
#print solver.solve(u, fixed_dofs)
#plot(u, interactive=True)
