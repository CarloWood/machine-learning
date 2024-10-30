#pragma once

#include "Direction.h"
#include "Point.h"

namespace math {

class Line
{
 protected:
  Point point_;
  Direction direction_;

 public:
  // Construct an undefined line.
  Line() = default;
  // Construct a line through point with direction.
  Line(Point const& point, Direction const& direction) : point_(point), direction_(direction) { }

  Point const& point() const { return point_; }
  Direction const& direction() const { return direction_; }

  operator Direction const&() const { return direction_; }

  Point intersection_with(Line const& line2) const;
};

} // namespace math
