from dolfin import *

mesh = Mesh("dolfin_coarse.xml")

outer_boundary = DomainBoundary()
facet_f = FacetFunction("size_t", mesh)
facet_f.set_all(0)
outer_boundary.mark(facet_f, 1) # gets the dolfin and the exterior

for facet in facets(mesh):
  if facet_f[facet]:
    x = facet.midpoint().x()
    y = facet.midpoint().y()
    if near(x, 0) or near(x, 1) or near(y, 0) or near(y, 1):
      facet_f[facet] = 0
plot(facet_f, interactive=True)


dolfin_vertices = []
for facet in facets(mesh):
  if facet_f[facet]:
    for vertex in vertices(facet):
      if vertex.index() not in dolfin_vertices:
        dolfin_vertices.append(vertex.index())
print dolfin_vertices

vertex_f = VertexFunction("size_t", mesh)
vertex_f.set_all(0)
for vertex in dolfin_vertices:
  vertex_f[vertex] = 1
plot(vertex_f, interactive=True)

mesh_coordinates = mesh.coordinates()



#for dolfin_vertex in dolfin_vertices:
#  print mesh_coordinates[]

