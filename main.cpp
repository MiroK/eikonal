#include "eikonal.h"
#include "test/CG1_VECTOR.h"
#include <dolfin.h>
#include <cstdlib>

using namespace eikonal;

// all available tests
template<typename T> int run_test(std::string type, std::size_t precision=1);

int main(int argc, char* argv[])
{
  dolfin::UnitSquareMesh mesh(10, 10);
  CG1::FunctionSpace V(mesh);
 
  dolfin::Function u(V);
  dolfin::Function du_dx(V);
  dolfin::Function du_dy(V);
  
  double _P[2] = {0.5, 0.5}; std::vector<double> P(_P, _P+2);
  MyPoint point(P);
  Problem problem(point);
  problem.exact_solution(u, du_dx, du_dy);

  dolfin::plot(u);
  dolfin::interactive(true);

  CG1_VECTOR::FunctionSpace VV(mesh);
  dolfin::Function du(VV);
  std::vector<double> du_values(VV.dim());

  boost::shared_ptr<dolfin::GenericVector> dx_vector = du_dx.vector();
  boost::shared_ptr<dolfin::GenericVector> dy_vector = du_dy.vector();
  for(std::size_t i = 0; i < du_values.size()/2; i++)
  {
    du_values[2*i] = (*dx_vector)[i];
    du_values[2*i+1] = (*dy_vector)[i];
  }
 
  std::cout << VV.dim() << std::endl;
  std::cout << V.dim() << std::endl;

  du.vector()->set_local(du_values);
  dolfin::plot(du);
  dolfin::interactive(true);

  return 0;
}

//-----------------------------------------------------------------------------

template<typename T> int run_test(std::string type, std::size_t precision=1)
{
  std::cout << type << std::endl;
  if(type == std::string("line"))
  {
    double _A[2] = {-2., -2.}; std::vector<double> A(_A, _A + 2);
    double _B[2] = {2., -2}; std::vector<double> B(_B, _B+2);
    Segment segment("segment", A, B);
    Problem problem(segment);
    
    int status;
    // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
    RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
    status = linear_2D_test<T>(problem, mesh_gen0, precision, false);

    // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
    RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
    status = linear_2D_test<T>(problem, mesh_gen1, precision, false);
    
    // convergence test on meshes by gmsh 0 .. 6
    GmshMeshGenerator mesh_gen2(1, 7, "rectangle");
    status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
    
    return status;
  }

  if(type == std::string("twocircle"))
  {
    double _c1[2] = {-1., 0.}; std::vector<double> c1(_c1, _c1+2);
    double _c2[2] = {sqrt(1.5), 0}; std::vector<double> c2(_c2, _c2+2);
    TwoCircles two_circles("twocircle", c1, 0.5, c2, 0.5);
    Problem problem(two_circles);
    
    int status;
    // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
    RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
    status = linear_2D_test<T>(problem, mesh_gen0, precision, false);

    // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
    RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
    status = linear_2D_test<T>(problem, mesh_gen1, precision, false);
    
    // convergence test on meshes by gmsh 0 .. 6
    GmshMeshGenerator mesh_gen2(1, 7, "rectangle");
    status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
    
    return status;
  }

  if(type == std::string("triangle"))
  {
    double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
    Polygon polygon("triangle", c, 1, 3);
    Problem problem(polygon);
    
    int status;
    // convergence test on meshes by gmsh 0 .. 6
    GmshMeshGenerator mesh_gen2(1, 7, "rectangle");
    status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
    
    return status;
  }

  if(type == std::string("zalesak"))
  {
    double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
    const double R = 1, W = 0.25, L = 1.5;
    Zalesak zalesak(c, R, W, L, 1000);
    Problem problem(zalesak);
    
    int status;
    // convergence test on meshes by gmsh 0 .. 6
    GmshMeshGenerator mesh_gen2(1, 7, "rectangle");
    status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
    
    return status;
  }

  if(type == std::string("dolphin"))
  {

    Dolphin dolphin;
    Problem problem(dolphin);
    
    int status;
    // convergence test on meshes by gmsh 0 .. 6
    GmshMeshGenerator mesh_gen2(1, 7, "rectangle");
    status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
    
    return status;
  }

  if(type == std::string("point"))
  {
    double _P[2] = {0., 0.}; std::vector<double> P(_P, _P+2);
    MyPoint point(P);
    Problem problem(point);
    
    int status;
    // convergence test on meshes by gmsh 0 .. 6
    GmshMeshGenerator mesh_gen2(1, 7, "rectangle");
    status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
    
    return status;
  }
  
  return 1;
}
