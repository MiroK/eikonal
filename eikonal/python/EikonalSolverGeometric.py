from EikonalSolverBase import *
from math import sin, cos, asin, acos, sqrt, pi
from numpy import dot

class EikonalSolverGeometric(EikonalSolverBase):
  '''Linear Eikonal solver'''
  def __init__(self, V):
    EikonalSolverBase.__init__(self, V)


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
