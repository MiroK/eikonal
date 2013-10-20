from EikonalSolverBase import *
from LineSampler import LineSampler
from scipy.optimize import fminbound, fmin_bfgs
import connectivity as conn
import numpy as np
from numpy.linalg import norm
from numpy import array

class Quad(EikonalSolverBase, LineSampler):
  '''Quadratic Eikonal solver'''
  def __init__(self, V, choice):
    EikonalSolverBase.__init__(self, V)
    
    # get also the other connectivities
    self.cell_to_edge = conn.build_cell_to_edge(V)
    self.edge_to_cell = conn.build_edge_to_cell(self.cell_to_edge)
    self.edge_to_dof =\
    conn.build_edge_to_dof(V, self.cell_to_edge)
    self.dof_to_edge = conn.build_dof_to_edge(self.edge_to_dof)
    self.choice = choice

  def order_dofs(self, dofs):
    '''Expects three colinear points and orders them according to their
    coordinates.'''
    if len(dofs) == 3:
      a = dofs[0]
      b = dofs[1]
      c = dofs[2]
      A = self.dof_to_coordinate[a]
      B = self.dof_to_coordinate[b]
      C = self.dof_to_coordinate[c]

      t = np.dot(C-A, B-A)/np.dot(B-A, B-A)
      if t > 1: # A B C, no need to swap
        return dofs
      elif t < 0: # C A B 
        return [c, a, b]
      else: # A C B
        return [a, c, b]

  def set_dofs_in_cell(self, dof, cell):
    '''Extract set dofs from cell.'''
    set_dofs = []
    for edge in set(self.cell_to_edge[cell]).difference(self.dof_to_edge[dof]):
      count = 0
      set_edge_dofs = []
      for edge_dof in self.edge_to_dof[edge]:
        if self.dof_status[edge_dof] and edge_dof != dof:
          count += 1
          set_edge_dofs.append(edge_dof)
        else:
          break

      if count == 3:
         set_dofs.append(self.order_dofs(set_edge_dofs))
   
    return set_dofs
    
  def local_solver(self, unset_U, u, _set_dofs, xtol):
    ''' Local solver on a triangle. Return new value for u in dof=unset_U. '''
    # set dofs can be a [[], []], i.e coming from multiple edges
    if self.choice == 0:
      return self.sample(unset_U, u, _set_dofs, xtol)
    
    retvals = []
    for set_dofs in _set_dofs: # get the minima from edge
      U = self.dof_to_coordinate[unset_U]
      A = self.dof_to_coordinate[set_dofs[0]]
      B = self.dof_to_coordinate[set_dofs[1]]
      C = self.dof_to_coordinate[set_dofs[2]]
     
      fp = u[unset_U]
      fA = u[set_dofs[0]]
      fB = u[set_dofs[1]]
      fC = u[set_dofs[2]]

      def value(t): # function to be minimized
        P = (1-t)*A + t*C
        d2 = sqrt(np.dot((U-P).flat, (U-P).flat))
        d1 = (1-t)*(1-2*t)*fA + 4*t*(1-t)*fB + t*(2*t-1)*fC
        return d1 + d2 

      # perform quadratic minimization
      output = fminbound(value, 0, 1, xtol=xtol, full_output=True)
      t_min = output[0]
      value_min = output[1]
      flag = output[2]
      
      f_min = min([value_min, fp])
      exact = norm(U - array([1., 0.]))
      
      print "old value", fp
      print "dof=", unset_U, "set_dofs", set_dofs, "value_min", value_min
      print "exact", exact, "error", abs(exact-value_min) 

      retvals.append(f_min)
    # return minima from all the edges
    return min(retvals)


  def solve(self, u, fixed_dofs, xtol):
    '''Global Eikonal solver.'''
    # we are expecting u from V
    u_vector = u.vector()
    if len(u_vector) != len(self.nodes_to_set):
      raise ValueError("u in not from the same function space as the solver.")

    # modify dof_status using fixed_dofs
    for dof in fixed_dofs:
      self.dof_status[dof] = True

    # remove dofs that are in fixed_dofs from nodes_to_set
    
    print "fixed", fixed_dofs
    
    for i in sorted(fixed_dofs)[::-1]:
      self.nodes_to_set.pop(i)

    print "nodes to set", self.nodes_to_set
    # we use four corners of mesh as reference points get them
    mesh = u.function_space().mesh()
    mesh_coordinates = mesh.coordinates()
    x_min = min(mesh_coordinates[:, 0])
    x_max = max(mesh_coordinates[:, 0])
    y_min = min(mesh_coordinates[:, 1])
    y_max = max(mesh_coordinates[:, 1])

    #create the sweep order
    sweep = []
    for x in [x_min, x_max]:
      for y in [y_min, y_max]:
        point = np.array([x, y])
        partial_sweep = self.sort_dofs(point, 2, False)
        sweep += partial_sweep
        sweep += partial_sweep[::-1]
    plot(u, interactive=True)    
    print sweep
    # create the reference function, fill it with some big
    v = Function(u.function_space())
    v.vector()[:] = max(u_vector.array())
    
    # start sweeping the mesh
    n_sweeps = 0
    off = len(sweep)/8
    x = 0
    CONTINUE = True
    while CONTINUE:
      for dof in sweep:
        print "global solver", dof
        x += 1
        u_old = u_vector[dof]
        for cell in self.dof_to_cell[dof]:
          _set_dofs = self.set_dofs_in_cell(dof, cell)
          print _set_dofs
          if _set_dofs:
            u_new =self.local_solver(dof, u_vector, _set_dofs, xtol)
            if u_new < u_old:
              self.dof_status[dof] = True
              print "dof status of", dof, self.dof_status[dof]
              u_old = u_new
        u_vector[dof] = u_old
        
        if x == off:
          n_sweeps += 1
          x = 0
          e = np.linalg.norm(u.vector().array() - v.vector().array(), np.inf)
          if e < DOLFIN_EPS:
            CONTINUE = False
            break
          else:
            v.vector()[:] = u.vector()[:]

    return n_sweeps
