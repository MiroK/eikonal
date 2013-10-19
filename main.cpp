#include "eikonal.h"
#include <cstdlib>

using namespace eikonal;

// all available tests
template<typename T> int run_test(std::string type, std::size_t precision=1);

int main(int argc, char* argv[])
{
  enum SOLVERS { LIN_GEOMETRIC, LIN_BRENT, LIN_NEWTON};

  if(argc == 3)
  {
    int solver = atoi(argv[1]);
    int precision = atoi(argv[2]);
    if(solver == LIN_GEOMETRIC)
    {
      std::cout << "Solving with linear geometric solver:" << std::endl;
      run_test<Solver>("line");
      run_test<Solver>("point");
      run_test<Solver>("two_circle");
      run_test<Solver>("triangle");
      run_test<Solver>("zalesak");
      run_test<Solver>("dolphin");
    }
    if(solver == LIN_BRENT)
    {
      std::cout << "Solving with linear Brent solver" << std::endl;
      std::cout << precision << std::endl;
      run_test<LinMinSolver>("line", precision);
      run_test<LinMinSolver>("point", precision);
      run_test<LinMinSolver>("two_circle", precision);
      run_test<LinMinSolver>("triangle", precision);
      run_test<LinMinSolver>("zalesak", precision);
      run_test<LinMinSolver>("dolphin", precision);
    }

    if(solver == LIN_NEWTON)
    {
      std::cout << "Solving with linear Newton solver:" << std::endl;
      std::cout << precision << std::endl;
      run_test<LinNewtonSolver>("line", precision);
      run_test<LinNewtonSolver>("point", precision);
      run_test<LinNewtonSolver>("two_circle", precision);
      run_test<LinNewtonSolver>("triangle", precision);
      run_test<LinNewtonSolver>("zalesak", precision);
      run_test<LinNewtonSolver>("dolphin", precision);
    }
  }
  return 1;
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

  if(type == std::string("two_circle"))
  {
    double _c1[2] = {-1., 0.}; std::vector<double> c1(_c1, _c1+2);
    double _c2[2] = {sqrt(1.5), 0}; std::vector<double> c2(_c2, _c2+2);
    TwoCircles two_circles("two_circles", c1, 0.5, c2, 0.5);
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
