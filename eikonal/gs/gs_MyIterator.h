#ifndef _GS_MY_ITERATOR_H_
#define _GS_MY_ITERATOR_H_

#include <vector>

namespace eikonal
{

  template<typename T>
  class MyIterator
  {
  /*
    Gives elements of std::vector<T> as if accessed by constant forward or
    reverse iterator without need to declare *::const_iterator or
    *::const_reverse_iterator.
  */
  public:
    // constructor, v=object that should be iterated, use reverse iterator
    MyIterator(const std::vector<T>& v, bool reverse=false) : start(v.begin()),
                                                              stop(v.end()),
                                                              current(start)
    {
      _next = not reverse ? &MyIterator::next : &MyIterator::reverse_next;
    }

    // copy-constructor, all variables are set according to other, especially
    // position of current
    MyIterator(const MyIterator& other) : start(other.start), stop(other.stop),
                                          current(other.current), 
                                          _next(other._next) { }

    // dereference
    T operator*() const { return ((this->*_next)()); }

    // increment
    void operator++() { current++; }

    // was the end of iteration reached
    bool end() const { return current == stop; }

  private:
    // next element in forward iteration
    T next() const { return *current; }

    // next element in reverse iteration
    T reverse_next() const
    { return *(stop - 1 - std::distance(start, current)); }

  private:
    // first element
    const typename std::vector<T>::const_iterator start;

    // last element
    const typename std::vector<T>::const_iterator stop;

    // current element;
    typename std::vector<T>::const_iterator current;

    // switch between forward and reverse iteration
    T (MyIterator::*_next)() const;
  };
}

#endif // _GS_MY_ITERATOR_H_
