#include "sys.h"
#include "Direction.h"
#include "Line.h"

namespace math {

Direction::Direction(Line const& line) : Direction(line.direction())
{
}

//static
Direction const Direction::up{0.0, 1.0};
//static
Direction const Direction::down{0.0, -1.0};
//static
Direction const Direction::left{-1.0, 0.0};
//static
Direction const Direction::right{1.0, 0.0};

#ifdef CWDEBUG
void Direction::print_on(std::ostream& os) const
{
  os << "{x_:" << x_ << ", y_:" << y_ << '}';
}
#endif

} // namespace math
