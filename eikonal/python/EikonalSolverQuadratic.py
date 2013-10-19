from EikonalSolverBase import *
from scipy.optimize import fminbound, fmin_bfgs
import connectivity as conn
import numpy as np
import matplotlib.pyplot as plt

def exact(U):
  x0_A = -1.
  y0_A = 0.

  x0_B = sqrt(1.5)
  y0_B = 0.

  r = 0.5
  
  a = abs(sqrt((U[0] - x0_A)**2+(U[1] - y0_A)**2) - r)
  b = abs(sqrt((U[0] - x0_B)**2+(U[1] - y0_B)**2) - r)
  return  min(a, b)


class EikonalSolverQuadratic(EikonalSolverBase):
  '''Quadratic Eikonal solver'''
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
    
  def local_solver(self, unset_U, u, _set_dofs, xtol):
    ''' Local solver on a triangle. Return new value for u in dof=unset_U. '''
    # set dofs can be a [[], []], i.e coming from multiple edges
    retvals = []
    #print "ls, unset_U,", unset_U, ", _set_dofs", _set_dofs
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

      plt.figure()
      t = np.linspace(-1, 1, 100)
      x = np.zeros(100)
      for i in range(len(x)):
        x[i] = value(t[i])
      plt.plot(t, x)
      plt.show()
     
      #plt.plot(t, x, label="quad")
      #plt.plot(t, y, label="lin")
      #plt.legend()
      #plt.show()
      # use the linear solver as initial guess
      #def value_linear(t):
      #  P = A*(1-t) + C*t
      #  return fA*(1-t) + fC*t + sqrt(np.dot((P-U).flat, (P-U).flat))
      #lin_out = fminbound(value_linear, 0, 1, full_output=True)
      #res = lin_out[0]
      #guess = lin_out[1]
      #lin_flag = lin_out[2]

      # perform quadratic minimization
      output = fminbound(value, 0, 1, xtol=xtol, full_output=True)
      t_min = output[0]
      value_min = output[1]
      flag = output[2]
      
      #print "linear", guess, "quadratic", value_min, "exact", exact(U)
      #print "using", fA, fB, fC, A, B, C
      
      # comparison with edge values
      a = sqrt(np.dot(U-A, U-A)) + fA
      b = sqrt(np.dot(U-C, U-C)) + fC


      #f_min = min([fA, fB, fC])
      f_min = min([value_min, fp])
      exact = sqrt((U[0] - 1)**2 + (U[1] - 0)**2)
      #print "fmi n", f_min, exact, "error=", abs(f_min - exact)
      retvals.append(f_min)
    #print "final", min(retvals)
    return min(retvals)
      #print "linear", guess, "quadratic", value_min, "exact", exact(U)
      #print "using", fA, fB, fC, A, B, C
      #print "crude", a, b 

      #if guess < exact(U) or value_min < exact(U):
      #  raise ValueError("NO!!!!!")

   
      #if flag == 0:
        #if value_min >= f_min:
        #  retvals.append(min([fp, value_min]))
        #else:
        #  if lin_flag == 0 and guess >= f_min: # use linear as fallback
        #    retvals.append(min([fp, guess]))
        #  else:
            #retvals.append(min([fp, a, b]))
            #print "dofs", A, B, C, "updating", U
            #print "using values", fA, fB, fC, fp
            #print "quadratic fo got", value_min
            #print "linear got", guess
            #print "derivative",  d_value(res), "hessian", dd_value(res)
            #res = fminbound(value, 0, 1)
            #print "linear got", res
            #t = np.linspace(0, 1, 100)
            #
            #plt.figure()      # plot nodes
            #plt.plot(A[0], A[1], "x" ,label="A")
            #plt.plot(B[0], B[1], "x" ,label="B")
            #plt.plot(C[0], C[1], "x" ,label="C")
            #plt.plot(U[0], U[1], "x" ,label="A")
            #plt.legend()
            #plt.show()

            #plt.figure()
            #t = np.linspace(-1, 1, 100)
            #x = np.zeros(100)
            #y = np.zeros(100)
            #for i in range(len(x)):
            #  x[i] = value(t[i])
            #  y[i] = value_linear(t[i])
           
            #plt.plot(t, x, label="quad")
            #plt.plot(t, y, label="lin")
            #plt.legend()
            #plt.show()
            #raise ValueError("Pushing value smaller than f_min") 
      #else:
      #    retvals.append(min([fp, a, b]))
    
    # return minimizer over all edges
