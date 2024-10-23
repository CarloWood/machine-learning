#pragma once

#include "utils/has_print_on.h"
#include <cmath>
#include <limits>
#include "debug.h"

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class CubicPolynomial;

class AnalyzedCubic
{
 protected:
  double inflection_point_;
  double signed_sqrt_D_{std::numeric_limits<double>::quiet_NaN()};
  double critical_point_w_;             // If extreme, then the one passed to initialize.
#if CW_DEBUG
  CubicPolynomial const* debug_cubic_{nullptr};
#endif

 public:
  void initialize(math::CubicPolynomial const& cubic, int extreme_type); // extreme_type: 1 (maximum) or -1 (minimum).

  bool has_extrema() const
  {
    // signed_sqrt_D_ is left as NaN when it is zero.
    return !std::isnan(signed_sqrt_D_);
  }

  double get_extreme() const
  {
    // Only call this function if has_extrema() returns true.
    ASSERT(has_extrema());
    return critical_point_w_;
  }

  double inflection_point() const
  {
    return inflection_point_;
  }

  double get_other_extreme() const
  {
    return 2.0 * inflection_point_ - critical_point_w_;
  }

#if CW_DEBUG
  // Return a pointer to the cubic that initalize was called with.
  math::CubicPolynomial const* debug_cubic() const { return debug_cubic_; }

  void print_on(std::ostream& os) const
  {
    os << "{inflection_point:" << inflection_point_ <<
         ", signed_sqrt_D:" << signed_sqrt_D_ <<
         ", critical_point_w:" << critical_point_w_ << '}';
  }
#endif
};

} // namespace math
