#pragma once

#include "LinePiece.h"
#include "Direction.h"
#include <cmath>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Vector
{
 protected:
  double x_;
  double y_;

 public:
  // Construct an uninitialized Vector.
  Vector() = default;

  // Construct a vector from its x,y coordinates.
  Vector(double x, double y) : x_(x), y_(y) { }

  // Construct a Vector that points in direction and has length.
  // Also used for automatic conversion from a Direction to a Vector.
  Vector(Direction direction, double length = 1.0) : x_(direction.x() * length), y_(direction.y() * length) { }

  // Construct a Vector from two points. If the second point is not given it defaults to the origin.
  // The direction is from the second argument (or origin) to the first argument.
  Vector(Point const& from, Point const& to) : x_(to.x() - from.x()), y_(to.y() - from.y()) { }
  explicit Vector(Point const& to) : x_(to.x()), y_(to.y()) { }

  // Construct a Vector from a LinePiece, pointing from the first point to the second point.
  Vector(LinePiece const& line_piece) : math::Vector(line_piece.from(), line_piece.to()) { }

  double x() const { return x_; }
  double y() const { return y_; }

  // Return dot product with v2.
  double dot(Vector const& v2) const { return x_ * v2.x_ + y_ * v2.y_; }

  // Return the cross product with v2.
  double cross(Vector const& v2) const { return x_ * v2.y_ - y_ * v2.x_; }

  // Construct a Direction from this vector.
  Direction direction() const { return Direction{Point{x_, y_}}; }

  // Return the length of the vector.
  double length() const { return std::hypot(x_, y_); }

  // Return the square of the length of the vector.
  double length_squared() const { return x_ * x_ + y_ * y_; }

  // Convert the vector to a Point.
  Point as_point() const { return {x_, y_}; }

 public:
  // Return the vector rotated 90 degrees counter-clockwise.
  Vector rotate_90_degrees() const { return { -y_, x_ }; }

  // Return the vector rotated 180 degrees.
  Vector rotate_180_degrees() const { return { -x_, -y_ }; }

  // Return the vector rotated 270 degrees.
  Vector rotate_270_degrees() const { return { y_, -x_ }; }

  // Add another vector.
  Vector& operator+=(Vector const& v2)
  {
    x_ += v2.x_;
    y_ += v2.y_;
    return *this;
  }

  // Subtract another vector.
  Vector& operator-=(Vector const& v2)
  {
    x_ -= v2.x_;
    y_ -= v2.y_;
    return *this;
  }

  // Multiply the vector with a scalar.
  Vector& operator*=(double scalar)
  {
    x_ *= scalar;
    y_ *= scalar;
    return *this;
  }

  // Divide the vector by a scalar.
  Vector& operator/=(double scalar)
  {
    x_ /= scalar;
    y_ /= scalar;
    return *this;
  }

  // Divide by a scalar.
  Vector operator/(double scalar) const
  {
    return {x_ / scalar, y_ / scalar};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << '[' << x_ << ", " << y_ << ']';
  }
#endif
};

inline Vector operator*(double length, Vector const& v2)
{
  return {length * v2.x(), length * v2.y()};
}

inline Point operator+(Point const& point, Vector const& v2)
{
  return {point.x() + v2.x(), point.y() + v2.y()};
}

inline Point operator-(Point const& point, Vector const& v2)
{
  return {point.x() - v2.x(), point.y() - v2.y()};
}

inline Vector operator+(Vector const& v1, Vector const& v2)
{
  return {v1.x() + v2.x(), v1.y() + v2.y()};
}

inline Vector operator-(Vector const& v1, Vector const& v2)
{
  return {v1.x() - v2.x(), v1.y() - v2.y()};
}

} // namespace math
