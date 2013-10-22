#ifndef __LINEAR_TEST_H
#define __LINEAR_TEST_H

/*
 all tests run on linear solvers
*/

#include "LinearTest.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

namespace eikonal
{
  class Segment;
  class Problem;
  class RectangleMeshGenerator;
  class GmshMeshGenerator;
  class Segment;
  class TwoCircles;
  class Polygon;
  class Zalesak;
  class MyPoint;
  class Dolphin;

  // type of object to compute signed distance from
  template<typename T>
  int run_linear_test(std::string type, std::size_t precision=1);
  
  // all tests
  int all_linear_tests(int solver, std::size_t precision);
}

// implementation

namespace eikonal
{
  template<typename T>
  int run_linear_test(std::string type, std::size_t precision=1)
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
      
      // convergence test on meshes by gmsh 0 .. 6 without smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 with smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = linear_2D_test<T>(problem, mesh_gen3, precision, false);
      
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
      
      // convergence test on meshes by gmsh 0 .. 6 without smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 with smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = linear_2D_test<T>(problem, mesh_gen3, precision, false);
      
      return status;
    }

    if(type == std::string("triangle"))
    {
      double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
      Polygon polygon("triangle", c, 1, 3);
      Problem problem(polygon);
      
      int status;
      // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
      status = linear_2D_test<T>(problem, mesh_gen0, precision, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = linear_2D_test<T>(problem, mesh_gen1, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 without smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 with smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = linear_2D_test<T>(problem, mesh_gen3, precision, false);
      
      return status;
    }

    if(type == std::string("zalesak"))
    {
      double _c[2] = {0., 0.}; std::vector<double> c(_c, _c + 2);
      const double R = 1, W = 0.25, L = 1.5;
      Zalesak zalesak(c, R, W, L, 1000);
      Problem problem(zalesak);
      
      int status;
      // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
      status = linear_2D_test<T>(problem, mesh_gen0, precision, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = linear_2D_test<T>(problem, mesh_gen1, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 without smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 with smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = linear_2D_test<T>(problem, mesh_gen3, precision, false);
      
      return status;
    }

    if(type == std::string("dolphin"))
    {

      Dolphin dolphin;
      Problem problem(dolphin);
      
      int status;
      // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
      status = linear_2D_test<T>(problem, mesh_gen0, precision, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = linear_2D_test<T>(problem, mesh_gen1, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 without smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 with smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = linear_2D_test<T>(problem, mesh_gen3, precision, false);
      
      return status;
    }

    if(type == std::string("point"))
    {
      double _P[2] = {0., 0.}; std::vector<double> P(_P, _P+2);
      MyPoint point(P);
      Problem problem(point);
      
      int status;
      // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
      status = linear_2D_test<T>(problem, mesh_gen0, precision, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = linear_2D_test<T>(problem, mesh_gen1, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 without smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = linear_2D_test<T>(problem, mesh_gen2, precision, false);
      
      // convergence test on meshes by gmsh 0 .. 6 with smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = linear_2D_test<T>(problem, mesh_gen3, precision, false);
      
      return status;
    }
    
    return 1;
  }
}

#endif // __LINEAR_TEST_H
