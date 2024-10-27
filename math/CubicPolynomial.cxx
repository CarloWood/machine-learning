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

int CubicPolynomial::get_roots(std::array<double, 3>& roots_out, int& iterations) const
{
  DoutEntering(dc::notice, "CubicPolynomial::get_roots() for " << *this);

  if (coefficients_[3] == 0.0)
  {
    // The cubic is actually a quadratic.
    QuadraticPolynomial qp(coefficients_[0], coefficients_[1], coefficients_[2]);
    return qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[0]));
  }

#if 0
  // Step one: divide all coefficients by coefficients_[3]. This does not change the roots.
  double c0 = coefficients_[0] / coefficients_[3];
  double c1 = coefficients_[1] / coefficients_[3];
  double c2 = coefficients_[2] / coefficients_[3];

  // Transform the cubic into
  //
  //   Q(u) = C0 + u (u^2 - 3)
  //
  double d = utils::square(c2) - 3.0 * c1;
  double C0 = (27.0 * c0 - (utils::square(c2) + 3 * d) * c2) / (d * std::sqrt(d));
#endif

  // Assuming the cubic has local extremes, we have the following eight possibilities:
  //
  //    m_.coefficient[3] > 0         |      m_.coefficient[3] < 0
  //                                  |
  // Case A:                          |   Case E:                           ⎫
  //      .-.     /                   |      \     .-.                      ⎮
  //     /   \   /                    |       \   /   \                     ⎮
  //    /     `-´                     |        `-´     \                    ⎮
  // --O-↑---↑-↑----------> x-axis    |   ------↑-↑---↑-O------> x-axis     ⎮   ← E is a minimum > 0
  //     x   I E   x = I - 2(E - I)   |         E I   x   x = I + 2(I - E)  ⎮
  //                                  |                                     ⎬ Iy > 0
  // Case B:                          |   Case F:                           ⎮
  //     .-.       /                  |     \       .-.                     ⎮
  //    /   \     /                   |      \     /   \                    ⎮
  // --O-↑---\---/-------> x-axis     |   ----\---/-----O------> x-axis     ⎮   ← E is a minimum < 0
  //  /  x  I `E´  x = same as case A |        `E´ I   x \x = same as case E⎭
  //                                  |
  // Case C:                          |   Case G:                           ⎫
  //     .-.       /                  |     \       .-.                     ⎮
  // ---/-↑-\----↓O-------> x-axis    |   ---O-----/-↑-\-------> x-axis     ⎮   ← E is a maximum > 0
  //   /  E I\   x x = same as case D |      x\   /I E  \ x = same as case H⎮
  //  /       `-´                     |        `·´       \                  ⎮
  //                                  |                                     ⎬ Iy < 0
  // Case D:        /                 |   Case H:                           ⎮
  // --------------O------> x-axis    |   --O------------------> x-axis     ⎮   ← E is a maximum < 0
  //      .-.     /                   |      \     .-.                      ⎮
  //     /   \   /                    |       \   /   \                     ⎮
  //    /  ↑ ↑`-´↑                    |       ↑`·´↑ ↑  \                    ⎮
  //       E I   x  x = I + 2(I - E)  |       x   I E     x = I - 2(E - I)  ⎭
  //
  // Where E is the extreme returned by acubic.get_extreme() and I is the inflection point.
  // O is the root that we will try to find using Halley's method, using a starting point of x.
  // The starting point is set on the slope that crosses the x-axis: passed the local extreme that is the closest to O.
  // Note: in all cases the initial starting point x = 3Ix - 2Ex.

  // Determine if the inflection point is above or below the x-axis.
  bool inflection_point_y_larger_than_zero =
    13.5 * coefficients_[0] + coefficients_[2] * (utils::square(coefficients_[2] / coefficients_[3]) - 4.5 * (coefficients_[1] / coefficients_[3])) > 0.0;

  AnalyzedCubic acubic;
  // Calculate the minimum (-1) if the inflection point lays above the x-axis, and the maximum otherwise.
  // This way acubic.get_extreme() (E above) becomes the extreme that is the closest to the x-axis.
  acubic.initialize(*this, inflection_point_y_larger_than_zero ? -1 : 1);

  // Obtain the calculated inflection point.
  double const inflection_point_x = acubic.inflection_point();

  // Remember if we have local extrema or not.
  bool const cubic_has_local_extrema = acubic.has_extrema();

  // Special case for if zero is a root (i.e. p(x) = x * (b + c x + d x^2)).
  double x = 0.0;

  if (AI_LIKELY(coefficients_[0] != 0.0))       // Is 0 not a root of the cubic?
  {
    // Avoid the local extrema and the inflection point because the derivative might be zero there too.
    x = cubic_has_local_extrema ?
      3 * inflection_point_x - 2 * acubic.get_extreme() :
      inflection_point_x + (((coefficients_[3] > 0.0) == inflection_point_y_larger_than_zero) ? -1.0 : 1.0);

    int limit = 100;
    double prev_x;
    double step = std::numeric_limits<double>::infinity();
    double prev_step;
    do
    {
      prev_x = x;
      prev_step = step;
      double f_x = this->operator()(x);                                         // c₀ + c₁x + c₂x² + c₃x³
      double three_c3_x = 3.0 * coefficients_[3] * x;                           // 3c₃x
      double half_fpp_x = coefficients_[2] + three_c3_x;                        // ½ ∂²f/∂x² = ½(2c₂ + 6c₃x) = c₂ + 3c₃x
      double fp_x = coefficients_[1] + (coefficients_[2] + half_fpp_x) * x;     // ∂f/∂x = c₁ + 2c₂x + 3c₃x² = c₁ + (2c₂ + 3c₃x)x
      // Apply Halley's method.
      step = -f_x * fp_x / (utils::square(fp_x) - f_x * half_fpp_x);
      x += step;                                                                // xₙ₊₁ = xₙ - f(x)f'(x) / (f'(x)² - ½f(x)f"(x))
      Dout(dc::notice, "Halley: x = " << std::setprecision(15) << x << "; Δx = " << step);
    }
    while (step != 0.0 && std::abs(step) < std::abs(prev_step) && --limit);
    iterations = 100 - limit;
  }

  roots_out[0] = x;
  int number_of_roots = 1;

  if (cubic_has_local_extrema)
  {
    // Find the other two roots, if any.
    [[maybe_unused]] double remainder;
    QuadraticPolynomial qp = long_division(roots_out[0], remainder);
    number_of_roots += qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[1]));
  }

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
