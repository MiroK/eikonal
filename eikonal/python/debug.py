# coding: utf-8
from Seeder import TwoCircle
from numpy import array
A = array([0, 0])
B = array([1, 0])
segment = TwoCircle(A, 0.5, B, 0.5)
from dolfin import *
mesh = RectangleMesh(-1., -1., 2., 1., 30, 30)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
from Problem import Problem
problem = Problem(segment)
problem.exact_solution(u)
plot(u, interactive=True)
