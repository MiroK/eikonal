#ifndef __HERMITE_TEST_H
#define __HERMITE_TEST_H

/*
 all tests run on hermite solvers
*/

#include "HermiteTest.h"
#include "Seeder.h"
#include "Problem.h"
#include "GmshMeshGenerator.h"
#include "RectangleMeshGenerator.h"
#include <string>

namespace eikonal
{
  // type of object to compute signed distance from
  template<typename T>
  int run_hermite_test(std::string solver_type,
                       std::size_t precision,
                       std::string ordering,
                       std::size_t p);
  
  // all tests
  int all_hermite_tests(std::size_t type, std::string ordering, std::size_t p);
}

// implementation
namespace eikonal
{
  template<typename T>
  int run_hermite_test(std::string type, std::size_t precision,
                       std::string ordering, std::size_t p)
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
      status = hermite_test<T>(problem, mesh_gen0, precision, ordering, p, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = hermite_test<T>(problem, mesh_gen1, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, no smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = hermite_test<T>(problem, mesh_gen2, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = hermite_test<T>(problem, mesh_gen3, precision, ordering, p, false);
      
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
      status = hermite_test<T>(problem, mesh_gen0, precision, ordering, p, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = hermite_test<T>(problem, mesh_gen1, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, no smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = hermite_test<T>(problem, mesh_gen2, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = hermite_test<T>(problem, mesh_gen3, precision, ordering, p, false);
      
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
      status = hermite_test<T>(problem, mesh_gen0, precision, ordering, p, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = hermite_test<T>(problem, mesh_gen1, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, no smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = hermite_test<T>(problem, mesh_gen2, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = hermite_test<T>(problem, mesh_gen3, precision, ordering, p, false);
      
      
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
      status = hermite_test<T>(problem, mesh_gen0, precision, ordering, p, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = hermite_test<T>(problem, mesh_gen1, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, no smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = hermite_test<T>(problem, mesh_gen2, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = hermite_test<T>(problem, mesh_gen3, precision, ordering, p, false);
      
      return status;
    }

    if(type == std::string("dolphin"))
    {

      Dolphin dolphin;
      Problem problem(dolphin);
      
      int status;
      // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen0(3, 8, -2, -2, 2, 2, false);
      status = hermite_test<T>(problem, mesh_gen0, precision, ordering, p, false);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(3, 8, -2, -2, 2, 2, true);
      status = hermite_test<T>(problem, mesh_gen1, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, no smoothing
      GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = hermite_test<T>(problem, mesh_gen2, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = hermite_test<T>(problem, mesh_gen3, precision, ordering, p, false);
      
      return status;
    }

    if(type == std::string("point"))
    {
      double _P[2] = {1., 0.}; std::vector<double> P(_P, _P+2);
      MyPoint point(P);
      Problem problem(point);
      
      int status;
      // convergence test on an unperturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen0(3, 7, 0, 0, 1, 1, false);
      status = hermite_test<T>(problem, mesh_gen0, precision, ordering, p, true);

      // convergence test on an perturbed RectangleMesh crossed [2**3 .. 2**7]
      RectangleMeshGenerator mesh_gen1(1, 3, 0, 0, 1, 1, true);
      status = hermite_test<T>(problem, mesh_gen1, precision, ordering, p, true);
      
      // convergence test on meshes by gmsh 0 .. 6, no smoothing
      /*GmshMeshGenerator mesh_gen2(1, 7, "rectangle", false);
      status = hermite_test<T>(problem, mesh_gen2, precision, ordering, p, false);
      
      // convergence test on meshes by gmsh 0 .. 6, smoothing
      GmshMeshGenerator mesh_gen3(1, 7, "rectangle", true);
      status = hermite_test<T>(problem, mesh_gen3, precision, ordering, p, false);
      */
      return status;
    }
    
    return 0;
  }
  //---------------------------------------------------------------------------
}
#endif // __HERMITE_TEST_H
