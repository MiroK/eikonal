'''
  This module includes functionality needed to set up the two circle test case.
  Solve the two-circle-problem (compute distance from circles (-1, 0, 0.5) and
  (sqrt(1.5), 0, 0.5) on [-2, 2]**2) using geometrical solver.
'''
from dolfin import *
import numpy

# problem parameters
x0_A = -1.
y0_A = 0.

x0_B = sqrt(1.5)
y0_B = 0.

r = 0.5
BIG_VALUE = 100

# exact solution
class ExactSolution(Expression):
  '''Distance from two circles.'''
  def eval(self, value, x):
    a = abs(sqrt((x[0] - x0_A)**2+(x[1] - y0_A)**2) - r)
    b = abs(sqrt((x[0] - x0_B)**2+(x[1] - y0_B)**2) - r)
    value[0] = min(a, b)

# initialize some function with correct values around 0 level set
# use markers for the job
class Particle:
  '''Lagrangian particle.'''
  def __init__(self, x, y):
    self.x = x
    self.y = y

  def __call__(self):
    return Point(self.x, self.y)

  def __str__(self):
    return "[%g, %g]" % (self.x, self.y)

  def position(self):
    return numpy.array([self.x, self.y])

def seed_circle(a, b, r, N):
  '''Seed circle (a, b, r) with N Lagrangian particles.'''
  theta = numpy.linspace(0, 2*numpy.pi, N)
  x = a + r*numpy.cos(theta)
  y = b + r*numpy.sin(theta)

  particles = []
  for i in range(N):
    particles.append(Particle(x[i], y[i]))

  return particles

def init_circle(u, x0, y0, r):
  '''Set values of u to distance from circle (x0, y0, r).'''
  V = u.function_space()
  mesh = V.mesh()
  dofmap = V.dofmap()
  u_vector = u.vector()

  num_cells = mesh.num_cells()
  num_particles = 10*num_cells        # just guess
  particles = seed_circle(x0, y0, r, num_particles)

  intersected_cells = []
  for particle in particles:
    cell_index = mesh.intersected_cell(particle())
    if cell_index != -1:
      intersected_cells.append(cell_index)

  intersected_cells = set(intersected_cells)
  
  set_dofs = []
  for cell in intersected_cells:
    dof_indices = dofmap.cell_dofs(cell)
    dof_positions = dofmap.tabulate_coordinates(Cell(mesh, cell))

    for dof_index, dof_position in zip(dof_indices, dof_positions):
      set_dofs.append(dof_index)
      x = dof_position[0]
      y = dof_position[1]

      d = abs(numpy.sqrt((x - x0)**2 + (y - y0)**2) - r)
      u_vector[dof_index] = d
  
  return set_dofs

def two_circle_problem(mesh, order, solver_class, plot_on=False):
  '''Solve the two circle problem with a solver. Compare solution in a norm.'''
  # initialize u 
  V = FunctionSpace(mesh, "CG", order)
  u = Function(V)
  
  # set with too far
  u.vector()[:] = BIG_VALUE

  # set intersected
  fixed_dofsA = init_circle(u, x0_A, y0_A, r)
  fixed_dofsB = init_circle(u, x0_B, y0_B, r)
  fixed_dofs = fixed_dofsA + fixed_dofsB
  fixed_dofs = set(fixed_dofs)

  #print fixed_dofs
  # initialize the solver and solve
  solver = solver_class(V)
  
  if plot_on:
    plot(u)
    interactive(True)
 
  num_iterations = solver.solve(u, fixed_dofs)

  if plot_on:
    plot(u)
    plot(interpolate(ExactSolution(), V))
    interactive(True)

  dim = V.dolfin_element().space_dimension()
  N = mesh.num_cells()
  out = File("results/u_%d_%d.pvd" % (dim, N))
  out << u

  errors = numpy.zeros(3)
  for norm_type in range(3):
    if norm_type == 0: # infty norm
      exact = interpolate(ExactSolution(), V)
      e = u.vector().array() - exact.vector().array();
      error = numpy.linalg.norm(e, numpy.inf)

    elif norm_type == 1: # l^1 norm
      exact = interpolate(ExactSolution(), V)
      e = u.vector().array() - exact.vector().array();
      error = numpy.linalg.norm(e, 1)/V.dim()
    
    elif norm_type == 2:
      error = errornorm(ExactSolution(), u)
  
    errors[norm_type] = error

  return mesh.hmin(), errors, num_iterations

