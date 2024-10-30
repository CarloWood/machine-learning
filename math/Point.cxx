#include "sys.h"
#include "Point.h"
#include "Direction.h"
#include "Vector.h"

namespace math {

Point Point::operator+(Direction const& direction)
{
  return {x_ + direction.x(), y_ + direction.y()};
}

Point Point::operator+(Vector const& v)
{
  return {x_ + v.x(), y_ + v.y()};
}

Point Point::operator-(Vector const& v)
{
  return {x_ - v.x(), y_ - v.y()};
}

Point& Point::operator+=(Direction const& direction)
{
  x_ += direction.x();
  y_ += direction.y();
  return *this;
}

Point& Point::operator+=(Vector const& v)
{
  x_ += v.x();
  y_ += v.y();
  return *this;
}

Point& Point::operator-=(Vector const& v)
{
  x_ -= v.x();
  y_ -= v.y();
  return *this;
}

Vector operator-(Point const& to, Point const& from)
{
  return {from, to};
}

bool operator!=(Point const& p1, Point const& p2)
{
  return p1.x_ != p2.x_ || p1.y_ != p2.y_;
}

#ifdef CWDEBUG
void Point::print_on(std::ostream& os) const
{
  os << "{x_:" << x_ << ", y_:" << y_ << '}';
}
#endif

} // namespace math
