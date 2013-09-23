#ifndef _TEST_H_
#define _TEST_H_

#include <vector>
#include <string>

namespace dolfin
{
  class Function;
  template<typename T> class MeshFunction;
}

/*
  Classes for defining boundary conditions (\Gamma) of different shapes. They
  are to be used for testing Eikonal solvers. For most of shapes there is an
  option to switch between initializing into distance function (all positive
  values) and signed distance function (positive and negative values). All 
  shapes are meant to work on UnitSquareMesh(es) or UnitCubeMesh(es).
*/

class GenericSeeder
{
public:
  // Initialize u into distance function (type == "d") or signed distance
  // function (type == "sd"). The degrees of freedom that are in cells inter-
  // sected by points of shape are returned into fixed_dofs. Mesh_f is used 
  // to define a domain where the Eikonal solver should operate later.
  virtual void initialize(dolfin::Function& u,
                          dolfin::MeshFunction<bool>& mesh_f,
                          std::vector<std::size_t>& fixed_dofs,
                          std::string type) = 0;
};

class LowerBoundarySeeder : public GenericSeeder
{
public:
  // constructor; k is a width of band above y=0 where mesh_f is set to true
  LowerBoundarySeeder(std::size_t k=0) : _k(k) { }

  // y = 0
  void initialize(dolfin::Function& u,
                  dolfin::MeshFunction<bool>& mesh_f,
                  std::vector<std::size_t>& fixed_dofs,
                  std::string type);
private:
  std::size_t _k;
};

#endif // _TEST_H_
