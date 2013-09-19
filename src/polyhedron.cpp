#include "polyhedron.h"

using namespace std;

std::vector<double>
barycenter(const std::vector<double>& vertices, const std::size_t n)
{
  std::size_t dim = vertices.size()/n;
  std::vector<double> _barycenter(dim);
  for(std::size_t i = 0; i < n; i++)
  {
    for(std::size_t j = 0; j < dim; j++)
    {
      _barycenter[j] += vertices[i*dim + j];
    }
  }

  for(std::size_t j = 0; j < dim; j++)
    _barycenter[j] /= n;

  return _barycenter;
}


