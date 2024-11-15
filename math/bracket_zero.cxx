#include "sys.h"
#include "bracket_zero.h"
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include "debug.h"

namespace math {

// Find a zero crossing between left and right and return the double representation
// that is closest to the middle of the region where f returns 0.
//
// Let `low` be the largest double representation such that f(low) < 0.
// Let `high` be the smallest double representation such that f(high) > 0.
//
// Then there are two possibilities,
//
// 1) f() returns 0 for one or more values of x:
//
//       f(x)
//        ^
// f_high→|          ──
//      0→+---- ───── ---> x
// f_low →|   ──
//       -+    ^     ^
//             |  ^  |
//            low | high
//              result
//
// In this case it is extremely likely that -f_low = f_high, just a single resolution step
// and the best result is 0.5 * (high + low).
//
// 2) f() is never zero:
//
//       f(x)
//        ^
// f_high→|     ──
//        |
//      0→+---------> x
//        |     ^result
//        |
// f_low →|   ──
//       -+    ^^
//             |high
//           low
//
// In which case abs(f_low) might be unequal to abs(f_high) and the best result is
// the x value for which abs(f(x)) is the smallest, where x is either low or high.
//
// Prerequisite: f(left) and f(right) are not allowed to be zero and must have different sign.
//               f must be monotonic between left and right.
//
double bracket_zero(double const left, double const right, std::function<double(double)> const& f)
{
  double f_left = f(left);
  double f_right = f(right);

  double const sign_of_f_right = std::copysign(1.0, f_right);

  // f_left and f_right must have opposite signs.
  if (f_left == 0 || f_right == 0 || std::copysign(1.0, f_left) == sign_of_f_right)
    throw std::runtime_error("math::bracket_zero: (left, right) does not bracket a zero.");

  double low = left;            // Find the largest value where f(x) < 0.
  double high = right;          // Find the smallest value where f(x) > 0.

  bool low_high_compare = low < high;
  while (std::nextafter(low, right) != high)
  {
    double mid = low + (high - low) / 2;
    double f_mid = f(mid);

    if (!std::isfinite(f_mid))
      throw std::runtime_error("Function returned infinite or NaN");

    if (f_mid == 0.0)
    {
      {
        double highz = mid;
        while (std::nextafter(low, right) != highz)
        {
          double midz = low + (highz - low) / 2;
          double f_mid = f(midz);

          if (!std::isfinite(f_mid))
            throw std::runtime_error("Function returned infinite or NaN");

          if (f_mid == 0.0 || std::copysign(1.0, f_mid) == sign_of_f_right)
            highz = midz;
          else
            low = midz;
        }
      }

      {
        double lowz = mid;
        while (std::nextafter(lowz, right) != high)
        {
          double midz = lowz + (high - lowz) / 2;
          double f_mid = f(midz);

          if (!std::isfinite(f_mid))
            throw std::runtime_error("Function returned infinite or NaN");

          if (f_mid == 0.0 || std::copysign(1.0, f_mid) != sign_of_f_right)
            lowz = midz;
          else
            high = midz;
        }
      }

      break;
    }
    else if (std::copysign(1.0, f_mid) == sign_of_f_right)
      high = mid;
    else
      low = mid;
  }

  double f_low = f(low);
  double f_high = f(high);

  if (std::nextafter(low, right) == high)
    return std::abs(f_low) < std::abs(f_high) ? low : high;

  low = std::nextafter(low, right);
  high = std::nextafter(high, left);

  // Calculate theoretical zero crossing using straight line approximation.
  // Note that both f_low and f_high are typically very close to zero and
  // f_low / (f_high - f_low) will typically be of the order of -0.5, so it
  // makes sense to calculate that first.
  return low - (f_low / (f_high - f_low)) * (high - low);
}

} // namespace math
