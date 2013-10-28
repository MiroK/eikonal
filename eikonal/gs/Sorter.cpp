#include "Sorter.h"
#include <dolfin/la/GenericVector.h>
#include <algorithm>

using namespace dolfin;

namespace eikonal
{
  Sorter::Sorter(const GenericVector& _values) : values(_values) { }
  //---------------------------------------------------------------------------

  void Sorter::sort(std::vector<dolfin::la_index>& dofs) const
  {
      std::sort(dofs.begin(), dofs.end(), *this);
  }
  //---------------------------------------------------------------------------
   
  bool Sorter::operator()(const dolfin::la_index& i,
                          const dolfin::la_index& j) const
  {

    return values[i] < values[j];
  }
  //---------------------------------------------------------------------------
}
