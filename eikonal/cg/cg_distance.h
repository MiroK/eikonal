#ifndef _CG_DISTANCE_H_
#define _CG_DISTANCE_H_

/*
  (Signed) distance functions between point and other objects.
*/

#include <vector>
#include <string>
#include <utility>

namespace eikonal
{
  // distance between 2 points in R^n
  double point_point(const std::vector<double>& P,
                     const std::vector<double>& Q);
  
  // distance between 2 points in R^n and unit vector (P-Q)/distnace //TODO
  void point_point_gradient(const std::vector<double>& P,
                            const std::vector<double>& Q,
                            std::vector<double>& gradient);
  
  // distance between point in R^2 and a line edge At+B(1-t), t\in[0, 1],
  // A, B in R^2
  double point_edge(const std::vector<double>& P,
                       const std::vector<double>& A,
                       const std::vector<double>& B);

  // gradient of the distance function from edge
  void point_edge_gradient(const std::vector<double>& P,
                           const std::vector<double>& A,
                           const std::vector<double>& B,
                           std::vector<double>& gradient);
  // (signed) distance between point and a line AB (as above but t\in R),
  // line and point make a half space so there is option for signed distance
  // pair.first is the (signed) distnace
  // pair.second indicates whether intersect was found within segment AB
  std::pair<double, bool> 
  point_line(const std::vector<double>& P, const std::vector<double>& A,
                    const std::vector<double>& B, const std::string type);
  
  // return also interesect
  std::pair<double, bool> 
  point_line(const std::vector<double>& P, const std::vector<double>& A,
             const std::vector<double>& B, const std::string type,
             std::vector<double>& I);

  // return gradient of distance function to line
  void point_line_gradient(const std::vector<double>& P,
                           const std::vector<double>& A,
                           const std::vector<double>& B,
                           std::vector<double>& gradient);

  // (signed) distance between point in R^2 and polygon given by flatt. vertices
  double point_polygon(const std::vector<double>& P,
                       const std::vector<double>& vertices,
                       const std::string type);
  
  // gradient of the distance function to polygon
  void point_polygon_gradient(const std::vector<double>& P,
                              const std::vector<double>& vertices,
                              std::vector<double>& _gradient);

  // (signed) distance between circle with center in R^2 and radius
  double point_circle(const std::vector<double>& P,
                      const std::vector<double>& center, const double radius,
                      const std::string type);
  
  // (signed) distance between point and plane given by ABC, all in R^3
  // pair.first is the distance
  // pair.second indicates whether intersect was found withing triangle ABC
  std::pair<double, bool> point_plane(const std::vector<double>& P,
                                      const std::vector<double>& A,
                                      const std::vector<double>& B,
                                      const std::vector<double>& C,
                                      const std::string type);
  
  // distance between point and a triangular face ABC, all in R^3
  double point_3face(const std::vector<double>& P,
                     const std::vector<double>& A,
                     const std::vector<double>& B,
                     const std::vector<double>& C);
  
  // (signed) distance between point and a tetrahedon ABCD, all in R^3
  double point_tet(const std::vector<double>& P,
                   const std::vector<double>& vertices,
                   const std::string type);
  
  // (signed) distance between point and a sphere with center and rad, R^3
  double point_sphere(const std::vector<double>& P,
                      const std::vector<double>& center, const double radius,
                      const std::string type);
}
#endif // _CG_DISTANCE_H_
