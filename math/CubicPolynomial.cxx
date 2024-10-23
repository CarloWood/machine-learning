#include "sys.h"
#include "CubicPolynomial.h"
#include "AnalyzedCubic.h"

namespace math {

void CubicPolynomial::initialize(double x0, double y0, double dxdy0, double x1, double y1, double dxdy1)
{
  double delta_x_inverse = 1.0 / (x0 - x1);

  // The theory of this approach is described here:
  // https://math.stackexchange.com/a/4926903/489074
  double d = (dxdy0 + dxdy1 - 2.0 * (y0 - y1) * delta_x_inverse) * (delta_x_inverse * delta_x_inverse);
  double c = 0.5 * ((dxdy0 - dxdy1) * delta_x_inverse - 3.0 * d * (x0 + x1));
  double b = (x0 * dxdy1 - x1 * dxdy0) * delta_x_inverse + 3.0 * d * x0 * x1;

  coefficients_[3] = d;
  coefficients_[2] = c;
  coefficients_[1] = b;
  coefficients_[0] = 0.0;
  coefficients_[0] = y0 - operator()(x0);
}

int CubicPolynomial::get_extrema(std::array<double, 2>& extrema_out, bool left_most_first) const
{
  DoutEntering(dc::notice, "CubicPolynomial::get_extrema(extrema_out, left_most_first = " << std::boolalpha << left_most_first << ")");

  // The cubic is:
  // coefficients_[0] + coefficients_[1] * x + coefficients_[2] * x^2 + coefficients_[3] * x^3.
  // The derivative is:
  // coefficients_[1] + 2 * coefficients_[2] * x + 3 * coefficients_[3] * x^2.

  if (coefficients_[3] == 0.0)
  {
    // If coefficients_[3] is zero, then the derivative has a single root when
    // coefficients_[1] + 2 * coefficients_[2] * x = 0 -->
    // x = -coefficients_[1] / (2 * coefficients_[2]).
    extrema_out[0] = -0.5 * coefficients_[1] / coefficients_[2];
    return std::isfinite(extrema_out[0]) ? 1 : 0;
  }

  // The determinant is (2 * coefficients_[2])^2 - 4 * coefficients_[1] * (3 * coefficients_[3]).
  double one_fourth_D = utils::square(coefficients_[2]) - 3.0 * coefficients_[1] * coefficients_[3];
  Dout(dc::notice, "D = " << (4.0 * one_fourth_D));

  if (one_fourth_D < 0.0)
  {
    // The inflection point is
    // -(2 * coefficients_[2]) / (2 * (3 * coefficients_[3])) =
    // -coefficients_[2] / (3 * coefficients_[3]).

    // Write the inflection point to index 0.
    extrema_out[0] = -coefficients_[2] / (3.0 * coefficients_[3]);
    return 0;
  }

  // Use a sqrt with the same sign as coefficients_[2];
  double const signed_half_sqrt_D = std::copysign(std::sqrt(one_fourth_D), coefficients_[2]);
  // The first root that is calculated should go here (see brute_force_index_first_root.cxx).
  int index_first_root = (left_most_first ? (coefficients_[3] > 0) : false) == (std::copysign(1.0, coefficients_[2]) > 0.0);

  // Calculate the root closest to zero.
  extrema_out[index_first_root] = -coefficients_[1] / (coefficients_[2] + signed_half_sqrt_D);

  if (AI_UNLIKELY(std::isnan(extrema_out[0])))
  {
    // This means we must have divided by zero, which means that both, coefficients_[1] as well as sqrtD, must be zero.
    // The latter means that coefficients_[0] is zero (coefficients_[2] was already checked not to be zero).
    // Therefore we have: f(x) = c x^2 with one root at x=0.
    extrema_out[0] = 0.0;
    return 1;
  }

  // Calculate the root further away from zero.
  extrema_out[1 - index_first_root] = -(coefficients_[2] + signed_half_sqrt_D) / (3.0 * coefficients_[3]);

  Dout(dc::notice, "extrema_out = " << std::setprecision(std::numeric_limits<double>::digits10) << extrema_out);

  // The smallest one must be in index 0 if left_most_first is true.
  ASSERT(one_fourth_D == 0.0 || !left_most_first || extrema_out[0] <= extrema_out[1]);
  // The minimum must be in index 0 if left_most_first is false.
  ASSERT(one_fourth_D == 0.0 || left_most_first || (extrema_out[0] < extrema_out[1] == coefficients_[3] < 0.0));

  return (one_fourth_D == 0.0) ? 1 : 2;
}

int CubicPolynomial::get_roots(std::array<double, 3>& roots_out) const
{
  double inflection_point_x = inflection_point();
  double inflection_point_y = this->operator()(inflection_point_x);

  //    m_.coefficient[3] > 0            m_.coefficient[3] < 0
  //
  // Case A:                          Case E:                           ⎫
  //      .-.     /                      \     .-.                      ⎮
  //     /   \   /                        \   /   \                     ⎮
  //    /     `E´   x = I - 2(E - I)       `E´     \  x = I + 2(I - E)  ⎮   ← minimum > 0
  // --O------------------> x-axis    --------------O------> x-axis     ⎮
  //                                                                    ⎬ IP > 0
  // Case B:                          Case F:                           ⎮
  //     .-.       /                    \       .-.                     ⎮
  //    /   \     / x = I - 2(E - I)     \     /   \  x = I + 2(I - E)  ⎮
  // --O-----\---/-------> x-axis     ----\---/-----O------> x-axis     ⎮
  //  /       `E´                          `E´       \                  ⎭   ← minimum < 0
  //
  // Case C:                          Case G:                           ⎫
  //     .E.       /                    \       .E.                     ⎮   ← maximum > 0
  // ---/---\-----O-------> x-axis    ---O-----/---\-------> x-axis     ⎮
  //   /     \   /  x = I + 2(I - E)      \   /     \ x = I - 2(E - I)  ⎮
  //  /       `-´                          `·´       \                  ⎮
  //                                                                    ⎬ IP < 0
  // Case D:        /                 Case H:                           ⎮
  // --------------O------> x-axis    --O------------------> x-axis     ⎮
  //      .E.     / x = I + 2(I - E)     \     .E.    x = I - 2(E - I)  ⎮   ← maximum < 0
  //     /   \   /                        \   /   \                     ⎮
  //    /     `-´                          `·´     \                    ⎭
  //
  // Where O is the root that we will try to find using Newton-Raphson.
  // Note: in all cases the initial Newton-Raphson starting point is 3I - 2E.

  AnalyzedCubic acubic;
  acubic.initialize(*this, inflection_point_y > 0.0 ? -1 : 1);

  if (!acubic.has_extrema())
  {
    // Avoid the inflection point because the derivative might be zero there.
    double initial_x = inflection_point_x + (coefficients_[3] > 0.0 == inflection_point_y > 0.0 ? -1.0 : 1.0);
//    roots_out[0] = newton_raphson(initial_x);
    // There is only one root if the cubic doesn't have extrema.
    return 1;
  }

//  roots_out[0] = newton_raphson(2 * inflection_point_x - 3 * acubic.get_extreme());
  int number_of_roots = 1;

  [[maybe_unused]] double remainder;
  QuadraticPolynomial qp = long_division(roots_out[0], remainder);
  number_of_roots += qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[1]));

  return number_of_roots;
}

#ifdef CWDEBUG
void CubicPolynomial::print_on(std::ostream& os) const
{
  bool first = true;
  int exponent = 0;
  for (double coefficient : coefficients_)
  {
    if (coefficient != 0.0)
    {
      if (first)
        os << coefficient;
      else if (coefficient > 0.0)
        os << " + " << coefficient;
      else
        os << " - " << -coefficient;
      if (exponent > 0)
      {
        os << " x";
        if (exponent > 1)
          os << '^' << exponent;
      }
      first = false;
    }
    ++exponent;
  }
}
#endif

} // namespace math
