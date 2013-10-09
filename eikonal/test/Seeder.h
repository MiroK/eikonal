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
    const std::vector<double>& A;
    const std::vector<double>& B;
    const std::size_t dim;
  };
  
  // seed circle (center, radius) with num_points
 /* void seed_circle(const std::vector<double>& center,
                    const double radius,
                    const std::size_t num_points,
                    std::vector<dolfin::Point>& points);

  // seed circles (center1, radius1) and (ceter2, radius2) with num_points each
  void seed_two_circles(const std::vector<double>& center1,
                        const double radius1,
                        const std::vector<double>& center2,
                        const double radius2,
                        const std::size_t num_points,
                        std::vector<dolfin::Point>& points);

  // seed polygon with num_vertices incribed into circle (center, radius) 
  // with num_points
  void seed_polygon(const std::vector<double>& center,
                    const double radius,
                    const std::size_t num_vertices
                    const std::size_t num_points,
                    std::vector<dolfin::Point>& points);*/
}

#endif // _SEEDER_H_
