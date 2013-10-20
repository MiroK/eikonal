from LineSampler import LineSampler
from EikonalSolverBase import *
from scipy.optimize import fminbound, fmin_bfgs
import connectivity as conn
import numpy as np
from numpy.linalg import norm
from numpy import array
from dolfin import Cell

class QuadraticRatic(EikonalSolverBase, LineSampler):
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
 
  def sort_dofs(self, point, norm_type, reverse):
    '''Order nodes_to_set(dofs) by their l^norm_type distance from point. '''
    cells = range(len(self.cell_to_dof)) 
    
    def foo(cell_index): # get the midpoint as array
      cell = Cell(self.mesh, cell_index)
      midpoint = cell.midpoint()
      return array([midpoint.x(), midpoint.y()])
  
    d = point.shape[0]
    key = lambda c_index : norm(point - foo(c_index))

    sorted_dofs = []
    cells.sort( key=key, reverse = reverse) # sort cells
    for cell in cells: # add all dofs of the cell that are to be set
      for dof in self.cell_to_dof[cell]:
        if dof in self.nodes_to_set:
          sorted_dofs.append(dof)

    return sorted_dofs
  
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
      
      #print "old value", fp
      #print "dof=", unset_U, "set_dofs", set_dofs, "value_min", value_min
      #print "exact", exact, "error", abs(exact-value_min) 

      retvals.append(f_min)
    # return minima from all the edges
    return min(retvals)
