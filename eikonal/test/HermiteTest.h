#ifndef _HERMITE_TEST_H_
#define _HERMITE_TEST_H_

/*
  Perform convergence test of the Eikonal Hermite solver on different problems.
*/

#include <string>

namespace dolfin
{
  class Mesh;
}

namespace eikonal
{
  class Problem;
  class MeshGenerator;

  // TODO
  int hermite_test(Problem& problem, MeshGenerator& mesh_gen,
                   std::size_t precision, std::string ordering,
                   std::size_t p, bool plot_on=false);

  // TODO
  int hermite_test(Problem& problem, const dolfin::Mesh& mesh,
                                          std::size_t precision,
                                          std::string ordering,
                                          std::size_t p,
                                          std::size_t& num_iters,
                                          std::size_t& min_calls,
                                          std::size_t& max_calls,
                                          double& l1_norm,
                                          double& l2_norm,
                                          double& band_l1_norm,
                                          double& band_l2_norm,
                                          double& coo_norm,
                                          double& time,
                                          std::string u_file_name,
                                          std::string exact_file_name,
                                          std::string error_file_name,
                                          bool plot_on=false);
}

#endif // _HERMITE_TEST_H
