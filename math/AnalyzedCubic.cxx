#include "sys.h"
#include "CubicPolynomial.h"
#include "AnalyzedCubic.h"

namespace math {

void AnalyzedCubic::initialize(math::CubicPolynomial const& cubic, int extreme_type)
{
  // Use -1 for 'minimum' and 1 for 'maximum'.
  ASSERT(extreme_type == -1 || extreme_type == 1);

#if CW_DEBUG
  // Remember for which cubic this function was called.
  debug_cubic_ = &cubic;
#endif

  double half_sqrt_D;
  if (AI_UNLIKELY(cubic[3] == 0.0))
  {
    // We don't really have an inflection point.
    inflection_point_ = std::numeric_limits<double>::infinity();

    // In this case the cubic is a parabola:
    //
    //   A(w) = b w^2 + c w + d
    //
    // with derivative:
    //
    //   A'(w) = 2b w + c
    //
    // The second derivative is a constant:
    //
    //   A''(w) = 2b
    //
    // Hence whether or not r is a minimum or maximum depends on the sign of b.
    // If b < 0 then it is a maximum.

    if ((extreme_type == 1) == (cubic[2] < 0.0))
    {
      // The derivative has one root at -c / 2b.
      critical_point_w_ = -0.5 * cubic[1] / cubic[2];
      // See below, assuming cubic[3] = 0.
      half_sqrt_D = std::abs(cubic[2]);
    }
    else
    {
      // Pretend the critical_point_w_ is at plus infinity.
      critical_point_w_ = std::numeric_limits<double>::infinity();
      // We should never be using this to calculate a height...
      half_sqrt_D = std::numeric_limits<double>::infinity();
    }
  }
  else
  {
    inflection_point_ = cubic.inflection_point();

    double one_fourth_D = utils::square(cubic[2]) - 3.0 * cubic[1] * cubic[3];
    if (one_fourth_D <= 0.0)
    {
      // If the determinant is zero, then the cubic has no local extremes.
      return;
    }

    // Use a sqrt with the same sign as cubic[2];
    half_sqrt_D = std::sqrt(one_fourth_D);
    double const half_Q = std::copysign(half_sqrt_D, cubic[2]);

    // The roots of the derivative are:
    // x_0 = -c / (b + 0.5 * Q);
    // x_1 = (-b - 0.5 * Q) / (3 * a);
    // where
    // f''(x_0) = Q
    // f''(x_1) = -Q
    // Therefore if Q is positive then x_0 is the minimum and x_1 is the maximum
    // and if Q is negative then x_0 is the maximum and x_1 is the minimum.
    //
    // Note: if cubic[2] is zero (or close to zero due to floating point round of errors)
    // then the sign of half_Q is not well defined, but its absolute value is usually still
    // significant. However, in that case, x_0 == -x_1 and changing the sign of half_Q has
    // no real influence because that is exactly where we swap formula as well.
    if ((extreme_type == 1) == (half_Q < 0.0))
    {
      // Calculate the root closest to zero.
      critical_point_w_ = -cubic[1] / (cubic[2] + half_Q);
    }
    else
    {
      // Calculate the root further away from zero.
      critical_point_w_ = -(cubic[2] + half_Q) / (3.0 * cubic[3]);
    }
  }
  // Don't ask about the minus sign.
  signed_sqrt_D_ = -2 * extreme_type * half_sqrt_D;
}

} // namespace math
