#ifndef _SORTER_H_
#define _SORTER_H_

#include <vector>
#include <dolfin/common/types.h>

namespace dolfin
{
  class GenericVector;
}

namespace eikonal
{
  class Sorter // dofs by their values 
  {
  public:
    // constructor
    Sorter(const dolfin::GenericVector& _values);

    // the sort
    void sort(std::vector<dolfin::la_index>& dofs) const;
   
    // functor for sorting
    bool operator()(const dolfin::la_index& i, const dolfin::la_index& j) const;

  private:
    const dolfin::GenericVector& values;
  };
}

#endif // _SORTER_H_
