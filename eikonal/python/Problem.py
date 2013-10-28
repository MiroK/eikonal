from numpy import zeros

class Problem:
  '''Problem defined by seeder.'''
  def __init__(self, seeder):
    '''Set seeder.'''
    self.seeder = seeder;

  def init(self, u, fixed_dofs, only=False):
    '''Initiliaze values of u at fixed_dofs.'''
    points = self.seeder.seed(100)
    i_cells = u.function_space().mesh().intersected_cells(points)
    
    # set everything to far
    far = 10
    if not only:
      u.vector()[:] += far

    dofmap = u.function_space().dofmap()
    mesh = u.function_space().mesh()
    dim = u.function_space().dim()
    all_coordinates = dofmap.tabulate_all_coordinates(mesh).reshape(dim ,2)
    for i_cell in i_cells:
      cell_dofs = dofmap.cell_dofs(i_cell)

      for cell_dof in cell_dofs:
        if cell_dof not in fixed_dofs: # append to fixed_dofs
          fixed_dofs.append(cell_dof)
          #print all_coordinates[cell_dof]
          u.vector()[cell_dof] = self.seeder.distance(all_coordinates[cell_dof]) 
    #print u.vector().array()

  def exact_solution(self, u):
    '''Compute distance from object defined by seeder.'''
    dofmap = u.function_space().dofmap()
    dim = u.function_space().dim()
    mesh = u.function_space().mesh()
    
    values = zeros(dim)
    for i, point in enumerate(dofmap.tabulate_all_coordinates(mesh).reshape(dim,
    2)):
      values[i] = self.seeder.distance(point)
    u.vector().set_local(values)
