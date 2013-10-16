import connectivity as conn
from dolfin import assemble, dx, sqrt, DOLFIN_EPS, Function, plot
import numpy as np

class EikonalSolverBase:
  '''Base class of Eikonal solvers.'''

  def __init__(self, V):
    '''V is a function space from u (the solution of Eikonal equation) lives.'''
    # get the basic connectivity information
    self.cell_to_dof = conn.build_cell_to_dof(V)
    self.dof_to_cell = conn.build_dof_to_cell(self.cell_to_dof)

    # map between dofs and their physical coordinates
    self.dof_to_coordinate = conn.build_dof_to_coordinate(self.dof_to_cell, V)

    # for the solution we will need dof_status now all False 
    self.dof_status = dict((dof, False) for dof in self.dof_to_cell)
    
    # and nodes_to_set, now it is all nodes
    self.nodes_to_set = range(V.dim())

  def sort_dofs(self, point, norm_type, reverse):
    '''Order nodes_to_set(dofs) by their l^norm_type distance from point. '''
    d = point.shape[0]
    key = lambda dof :\
    sum([(abs(point[i] - self.dof_to_coordinate[dof][i]))**norm_type\
                            for i in range(d)])

    return sorted(self.nodes_to_set, key = key, reverse = reverse)

  def set_dofs_in_cell(self, dof, cell):
    '''
      Extract set dofs from cell. These are used for computing value in dof.
      Each solver needs different dofs, so this is implemented in children.
    '''
    raise NotImplementedError('''Implement in child!''')

  def local_solver(self,unset_U, u, set_dofs):
    ''' Local solver on a triangle. Return new value for u in dof=unset_U. '''
    raise NotImplementedError('''Implement in child!''')

  def solve(self, u, fixed_dofs):
    '''Global Eikonal solver.'''
    # we are expecting u from V
    u_vector = u.vector()
    if len(u_vector) != len(self.nodes_to_set):
      raise ValueError("u in not from the same function space as the solver.")

    # modify dof_status using fixed_dofs
    for dof in fixed_dofs:
      self.dof_status[dof] = True

    # remove dofs that are in fixed_dofs from nodes_to_set
    for i in sorted(fixed_dofs)[::-1]:
      self.nodes_to_set.pop(i)

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
    
    # create the reference function, fill it with some big
    v = Function(u.function_space())
    v.vector()[:] = max(u_vector.array())
    #print sweep    
    # start sweeping the mesh
    count = 0
    off = len(sweep)/8
    x = 0
    CONTINUE = True
    while CONTINUE:
      count += 1
      i = -off
      change_count = 0
      order = 0
      for dof in sweep:
        x += 1
        i += 1
        u_old = u_vector[dof]
        for cell in self.dof_to_cell[dof]:
          _set_dofs = self.set_dofs_in_cell(dof, cell)
          #print "Cell", cell, "dof", dof, _set_dofs, "u_old", u_old
          if _set_dofs:
            print "cell=", cell
            u_new =self.local_solver(dof, u_vector, _set_dofs)
            if u_new < 0:
              raise ValueError("negative value")
           # print "\t new value", u_new
            if u_new < u_old:
              change_count += 1
              self.dof_status[dof] = True
              u_old = u_new
              #print "u_old \t\t changed to", u_old 
        u_vector[dof] = u_old
        print "dof is",dof, "error= %.16f" % abs(u_old - self.dof_to_coordinate[dof][1]),\
        "num=%.18f" % u_old, "exact=", self.dof_to_coordinate[dof][1]
       # print "Final value", u_old
      
        if x == off:
          x = 0
          e = np.linalg.norm(u.vector().array() - v.vector().array(), np.inf)
          print x
          if e < DOLFIN_EPS:
            CONTINUE = False
          else:
            v.vector()[:] = u.vector()[:]
            order += 1
            i = -off


    return count
