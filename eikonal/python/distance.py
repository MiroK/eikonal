from numpy.linalg import norm
from numpy import dot

def point_line(P, A, B):
  '''Distance between point P and line defined by A, B.'''

  AB = norm(B-A)
  t = dot(P-A, B-A)/AB/AB

  inside = False if (t>1 or t<0) else True
  I = A + t*(B-A)

  return norm(P-I), inside
#-------------------------------------------------------------------------------

def point_edge(P, A, B):
  '''Distance between point P and edge A, B.'''

  (distance, inside) = point_line(P, A, B)
  return distance if inside else - 1
#-------------------------------------------------------------------------------

def point_circle(P, c, r):
  '''Distance between point P and circle(c, r). '''
  return abs(norm(P-c) - r)
#-------------------------------------------------------------------------------
