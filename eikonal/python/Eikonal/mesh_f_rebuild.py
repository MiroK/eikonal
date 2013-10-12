from dolfin import *

mesh = UnitSquareMesh(50, 50)

mesh_f = MeshFunction("bool", mesh, 2)
mesh_f.set_all(False)
for cell in cells(mesh):
  if 0.35 < cell.midpoint().y() < 0.65:
    mesh_f[cell] = True


V = FunctionSpace(mesh, "CG", 1)
dofmap = V.dofmap()
u = Function(V)
for cell in cells(mesh):
  if mesh_f[cell]:
    dofs = dofmap.cell_dofs(cell.index())
    for dof in dofs:
      u.vector()[dof] = 100

plot(u, interactive=True)

