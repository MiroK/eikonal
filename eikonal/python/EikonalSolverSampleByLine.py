from EikonalSolverBase import *
from scipy.optimize import fminbound, fmin_bfgs
import connectivity as conn
from math import sin, cos, asin, acos, sqrt, pi
from numpy import dot
import matplotlib.pyplot as plt


class EikonalSolverSampleByLine(EikonalSolverBase):
  def __init__(self, V):
    EikonalSolverBase.__init__(self, V)
    
    # get also the other connectivities
    self.cell_to_edge = conn.build_cell_to_edge(V)
    self.edge_to_cell = conn.build_edge_to_cell(self.cell_to_edge)
    self.edge_to_dof =\
    conn.build_edge_to_dof(V, self.cell_to_edge)
    self.dof_to_edge = conn.build_dof_to_edge(self.edge_to_dof)

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
    # it is rare bot possible that dof can be updated by two edges
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

  def lin_local_solver(self, unset_U, u, set_dofs, xtol):
    ''' Local solver on a triangle. Return new value for u in dof=unset_U. '''
    C = self.dof_to_coordinate[unset_U]
    A = self.dof_to_coordinate[set_dofs[0]]
    B = self.dof_to_coordinate[set_dofs[1]]
    
    u_A = u[set_dofs[0]]
    u_B = u[set_dofs[1]]
    u_C = u[unset_U]
    
    c = sqrt(dot(B-A, B-A))
    b = sqrt(dot(C-A, C-A))
    a = sqrt(dot(C-B, C-B))

    alpha = acos(dot(C-B, A-B)/a/c)
    beta = acos(dot(C-A, B-A)/b/c)
    
    if abs(u_B-u_A) <= c:
      theta = asin((u_B-u_A)/c)

      if (max(0, alpha-pi/2) <= theta <= pi/2-beta) or \
         (alpha-pi/2 <= theta <= min(0, pi/2-beta)):
         
         h = a*sin(alpha-theta)
         H = b*sin(beta+theta)

         u_ = 0.5*(h+u_B) + 0.5*(H+u_A)
         return min(u_C, u_)

      else:
        u_ = u_A + b
        _u = u_B + a
        return min([u_C, u_, _u])
    
    else:
      u_ = u_A + b
      _u = u_B + a
      return min([u_C, u_, _u])

  def local_solver(self, unset_U, u, _set_dofs, xtol):
    ''' Local solver on a triangle. Return new value for u in dof=unset_U. '''
    _set_dofs
    # set dofs can be a [[], []], i.e coming from multiple edges
    retvals = []
    for set_dofs in _set_dofs: # get the minima from edge
      one = self.lin_local_solver(unset_U, u, [set_dofs[0], set_dofs[1]], xtol)    
      two = self.lin_local_solver(unset_U, u, [set_dofs[1], set_dofs[2]], xtol)    
      tre = self.lin_local_solver(unset_U, u, [set_dofs[2], set_dofs[0]], xtol)    
      retvals.append(min([one, two, tre]))

    #print "final", min(retvals)
    return min(retvals)
