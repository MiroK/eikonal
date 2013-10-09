#include "eikonal.h"

using namespace eikonal;

int main()
{
  // define the seeder
  double _A[2] = {0, 0}; std::vector<double> A(_A, _A + 2);
  double _B[2] = {1, 0}; std::vector<double> B(_B, _B + 2);
  Segment segment("segment_[0,0]_[1,0]", A, B);

  // define the problem
  Problem problem(segment);
 
  // perform convergence test using Solver
  int status = linear_2D_test<Solver>(problem, 3, 5, false);

  return 0;
}
