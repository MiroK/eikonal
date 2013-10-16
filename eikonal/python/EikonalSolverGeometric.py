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

  def local_solver(self, unset_U, u, set_dofs):
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
    
    print "local solver, set_dofs=", set_dofs, "alpha, beta", alpha, beta
    print A, B, u_A, u_B

    if abs(u_B-u_A) <= c:
      theta = asin(abs(u_B-u_A)/c)
      print "c=", c, "theta=", theta
      cond1 = max(0, alpha-pi/2) <= theta <= pi/2-beta
      cond2 = alpha-pi/2 <= theta <= min(0, pi/2-beta)
      print cond1, cond2
      print max(0, alpha-pi/2), theta, pi/2 - beta
      if (max(0, alpha-pi/2) <= theta <= pi/2-beta) or \
         (alpha-pi/2 <= theta <= min(0, pi/2-beta)):
         
         theta = theta if u_A <= u_B else -theta;

         h = a*sin(alpha-theta)
         H = b*sin(beta+theta)

         u_ = 0.5*(h+u_B) + 0.5*(H+u_A)
         return min(u_C, u_)

      else:
        u_ = u_A + b
        _u = u_B + a
        print "Returning on alfa, beta fail."
        print "u_A=%.12f u_B=%.12f a=%.12f b=%.12f" %  (u_A, u_B, a, b)
        print "A=[%.12f, %.12f], B=[%.12f, %.12f] C[1]=%.12f" % \
        (A[0], A[1], B[0], B[1], C[1])
        print min([u_C, u_, _u])
        return min([u_C, u_, _u])
    
    else:
      u_ = u_A + b
      _u = u_B + a
      print "Returning on theta fail"
      return min([u_C, u_, _u])

