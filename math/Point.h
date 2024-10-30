#pragma once

#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Direction;
class Vector;

class Point
{
 protected:
  double x_;
  double y_;

 public:
  Point() = default;
  Point(double x, double y) : x_(x), y_(y) { }

  double x() const { return x_; }
  double y() const { return y_; }

  Point operator+(Direction const& direction);
  Point operator+(Vector const& v);
  Point operator-(Vector const& v);
  Point& operator+=(Direction const& direction);
  Point& operator+=(Vector const& v);
  Point& operator-=(Vector const& v);
  friend Vector operator-(Point const& to, Point const& from);
  friend bool operator!=(Point const& p1, Point const& p2);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace math
