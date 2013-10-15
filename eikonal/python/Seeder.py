from dolfin import Point
from distance import point_edge, point_circle
from numpy import linspace
from math import pi, cos, sin

class TwoCircle:
  '''Seeder for two circle problem.'''
  def __init__(self, c1, r1, c2, r2):
    '''Set properties of circle.'''
    self.c1 = c1
    self.r1 = r1
    self.c2 = c2
    self.r2 = r2

  def seed(self, num_points = 100):
    '''Put num_points points on each circle.'''
    thetas = linspace(0, 2*pi, num_points)

    points = []
    for theta in thetas:
      x = self.c1[0] + self.r1*cos(theta)
      y = self.c1[1] + self.r1*sin(theta)
      points.append(Point(x, y))
    
    for theta in thetas:
      x = self.c2[0] + self.r2*cos(theta)
      y = self.c2[1] + self.r2*sin(theta)
      points.append(Point(x, y))

    return points

  def distance(self, point):
    '''Compute distance of point from two circles.'''
    return min(point_circle(point, self.c1, self.r1),\
               point_circle(point, self.c2, self.r2))
#-------------------------------------------------------------------------------

class Segment:
  '''Seeder for segment.'''
  
  def __init__(self, A, B):
    '''Set the points defining segment.'''
    self.A = A
    self.B = B

  def seed(self, num_points=1000):
    '''Put num_points points on the segment.'''
    ts = linspace(0, 1, num_points)
    
    points = []
    for t in ts:
      P = self.A + t*(self.B-self.A)
      points.append(Point(float(P[0]), float(P[1])))

    return points

  def distance(self, P):
    '''Compute distance of point from segment.'''
    return point_edge(P, self.A, self. B)
#-------------------------------------------------------------------------------
