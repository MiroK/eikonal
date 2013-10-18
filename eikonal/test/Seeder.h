#ifndef _SEEDER_H_
#define _SEEDER_H_

/*
  Seed surface of an object with dolfin::Points
*/

#include <vector>
#include <string>

namespace dolfin
{
  class Point;
} 

namespace eikonal
{
  class Seeder // Interface fo interacting with Problem
  {
  public:
    // constructor, set name
    Seeder(const std::string& _name) : name(_name) { }

    // seed object with num_points points
    virtual void seed(std::vector<dolfin::Point>& points,
                      const std::size_t num_points) const = 0;

    // compute distance of point from object
    virtual double distance(const std::vector<double>& point) const = 0;
  
  public:
    const std::string name;
  };
  //---------------------------------------------------------------------------

  class Segment : public Seeder
  {
  public:
    // constructor, segment between point _A, _B
    Segment(const std::string& name,
            const std::vector<double>& _A,
            const std::vector<double>& _B);
    
    // put num_points points on the segment
    void seed(std::vector<dolfin::Point>& points,
              const std::size_t num_points) const;
    
    // compute distance from the segment 
    double distance(const std::vector<double>& point) const;  
    
  private:
    const std::vector<double> A;
    const std::vector<double> B;
    const std::size_t dim;
  };
  //----------------------------------------------------------------------------

  class TwoCircles : public Seeder
  {
  public:
    // constuctor, two circles with center c_i and radius r_i
    TwoCircles(const std::string& name,
              const std::vector<double>& _c1, const double _r1,
              const std::vector<double>& _c2, const double _r2);

    // put num_points on each circle
    void seed(std::vector<dolfin::Point>& points,
              const std::size_t num_points) const;
      
    // compute distance from the sigment, distance NOT SIGNED DISTANCE
    double distance(const std::vector<double>& point) const;

  private:
    const std::vector<double> c1;
    const std::vector<double> c2;
    const double r1;
    const double r2;
  };
  //----------------------------------------------------------------------------
 
  class Polygon : public Seeder
  {
  public:
    // constructor, polygon with n vertices, incribed into circle centered
    // at c with radius r
    Polygon(const std::string& name, const std::vector<double>& c,
            const double r, const std::size_t n);

    // put num_points on each edge of the polygon
    void seed(std::vector<dolfin::Point>& points,
              const std::size_t num_points) const;

    // compute distance from the polygon
    double distance(const std::vector<double>& point) const;

  private:
    std::vector<double> vertices;
  };
  //----------------------------------------------------------------------------

  class Zalesak : public Seeder
  {
  public:
    // constructor, Zalesak disk defined as polygon with by num_vertices
    // center, radius, notch length and notch width
    Zalesak(const std::vector<double>& c,
            const double R, const double W, const double L, 
            const std::size_t num_vertices);

    // here num_points are ignored and 3 points are put on each edge
    void seed(std::vector<dolfin::Point>& points,
              const std::size_t num_points) const;

    // compute distance from the disk
    double distance(const std::vector<double>& point) const;

  private:
    std::vector<double> vertices;
  };
  //----------------------------------------------------------------------------

  class Dolphin : public Seeder
  {
  public:
    // seed dolfin usin points in file 
    Dolphin();

    // here num_points are ignored and
    void seed(std::vector<dolfin::Point>& points,
              const std::size_t num_points) const;

    // compute distance from dolfin
    double distance(const std::vector<double>& point) const;

  private:
    std::vector<double> vertices;
  };
  //----------------------------------------------------------------------------
}
#endif // _SEEDER_H_
