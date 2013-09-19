#ifndef _POLYHEDRON_H_
#define _POLYHEDRON_H_

/*
  Functions for computations with with polyhedrons.
*/

#include <vector>

// compute barycenter of poly(gon/hedron) specified by vertices
// vertices are flattened so number of vertices must be given to avoid ambiguity
std::vector<double> barycenter(const std::vector<double>& vertices,
                               const std::size_t n);

#endif // _POLYHEDRON_H_
