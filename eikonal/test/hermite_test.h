#ifndef __HERMITE_TEST_H
#define __HERMITE_TEST_H

/*
 all tests run on hermite solvers
*/

#include "HermiteTest.h"
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
  int run_hermite_test(std::string type, std::size_t precision=1);
  
  // all tests
  int all_hermite_tests(std::size_t precision);
}


#endif // __HERMITE_TEST_H
