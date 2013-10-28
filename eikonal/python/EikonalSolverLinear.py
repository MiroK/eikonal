from EikonalSolverBase import *
from scipy.optimize import fminbound, fmin_bfgs

class EikonalSolverLinear(EikonalSolverBase):
  '''Linear Eikonal solver'''
  def __init__(self, V):
    EikonalSolverBase.__init__(self, V)
    self.max_calls = 0   # maximum number of calls to minimized function 
    self.min_calls = 501 # minimum ...

  def set_dofs_in_cell(self, dof, cell):
    ''' Extract set dofs from cell.'''
    set_dofs = []
    for _dof in self.cell_to_dof[cell]:
      if _dof != dof and self.dof_status[_dof]:
        set_dofs.append(_dof)

    if len(set_dofs) == 2:
      return set_dofs
    else:
      return []

  def local_solver(self, unset_U, u, set_dofs, xtol):
    ''' Local solver on a triangle. Return new value for u in dof=unset_U. '''
    p = self.dof_to_coordinate[unset_U]
    p0 = self.dof_to_coordinate[set_dofs[0]]
    p1 = self.dof_to_coordinate[set_dofs[1]]
    
    fp = u[unset_U]
    f0 = u[set_dofs[0]]
    f1 = u[set_dofs[1]]

    def value(t):
      P = p0*(1-t) + p1*t
      return f0*(1-t) + f1*t + sqrt(np.dot((P-p).flat, (P-p).flat))

    res, fval, ierr, numfunc = fminbound(value, 0, 1, xtol=xtol, full_output=True)
    
    # set the "records" of eval calls
    if numfunc > self.max_calls:
      self.max_calls = numfunc

    if numfunc < self.min_calls:
      self.min_calls = numfunc
    
    if ierr == 0:
      return min([value(res), fp])
    else:
      # edge lentgh
      a = sqrt(np.dot(p-p0, p-p0)) + f0
      b = sqrt(np.dot(p-p1, p-p1)) + f1

      return min([fp, a, b])
