#include "cg_common.h"

using namespace std;

namespace eikonal
{
  std::vector<double> barycenter(const std::vector<double>& vertices,
                                 const std::size_t dim)
  { 
    const size_t n_vertices = vertices.size()/dim;
    vector<double> _barycenter(dim);
    for(size_t x = 0; x < dim; x++)
    {
      for(size_t i = 0; i < n_vertices; i++)
      {
        _barycenter[x] += vertices[i*dim + x];
      }
      _barycenter[x] /= n_vertices;
    }
    return _barycenter;
  }
  //---------------------------------------------------------------------------

}
