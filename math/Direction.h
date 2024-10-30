#pragma once

#include "LinePiece.h"
#include <cmath>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif
#include "debug.h"

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Line;

class Direction
{
 protected:
  double x_;    // x,y is a unit vector pointing in the intended direction.
  double y_;

 public:
  // Construct an undefined Direction.
  Direction() = default;

  // Construct a Direction that points in the direction theta (in radians): an angle with the positive x-axis.
  Direction(double theta) : x_(std::cos(theta)), y_(std::sin(theta)) { }

  // Construct a Direction from two points. If the second point is not given it defaults to the origin.
  explicit Direction(Point const& from, Point const& to) : x_(to.x() - from.x()), y_(to.y() - from.y())
  {
    double len = std::sqrt(x_ * x_ + y_ * y_);
    x_ /= len;
    y_ /= len;
    ASSERT(!std::isnan(x_) && !std::isnan(y_));
  }

  // If only one point is give, the direction is from the origin to that point.
  explicit Direction(Point const& to) : Direction(Point{0.0, 0.0}, to) { }

  // Construct a Direction from a LinePiece, pointing from the first point to the second point.
  Direction(LinePiece const& line_piece) : Direction(line_piece.from(), line_piece.to()) { }

  // Construct a Direction from a Line.
  Direction(Line const& line);

  double x() const { return x_; }
  double y() const { return y_; }

  // Return dot product with d2.
  double dot(Direction const& d2) const { return x_ * d2.x_ + y_ * d2.y_; }

  // Returns an angle in the range (-π, π] radians.
  double as_angle() const { return std::atan2(y_, x_); }

 protected:
  // For normal() and inverse().
  constexpr Direction(double x, double y) : x_(x), y_(y) { };

 public:
  // Return the direction rotated 90 degrees counter-clockwise.
  Direction normal() const { return { -y_, x_ }; }

  // Return the direction rotated 180 degrees.
  Direction inverse() const { return { -x_, -y_ }; }

  // Return the direction rotated 270 degrees.
  Direction normal_inverse() const { return { y_, -x_ }; }

  static Direction const up;
  static Direction const down;
  static Direction const left;
  static Direction const right;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace math
