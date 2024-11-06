#include "sys.h"
#include "CubicPolynomial.h"
#include "AnalyzedCubic.h"
#ifdef CWDEBUG
#include "cwds/Restart.h"
#endif

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

  // Step one: divide all coefficients by coefficients_[3]. This does not change the roots.
  double const c0 = coefficients_[0] / coefficients_[3];
  double const c1 = coefficients_[1] / coefficients_[3];
  double const c2 = coefficients_[2] / coefficients_[3];
  double const d = utils::square(c2) - 3.0 * c1;

  // The cubic is now monic:
  //
  //   p(x) = c0 + c1 x + c2 x² + x³
  //
  // The first derivative is,
  //
  //   p'(x) = c1 + 2 c2 x + 3 x²
  //
  // Setting this to zero gives:
  //                                  -c2 +/- sqrt(c2^2 - 3c1)   -c2 +/- sqrt(d)
  //   0 = c1 + 2 c2 x + 3 x² --> x = ------------------------ = ---------------
  //                                             3                     3
  // Remember if we have local extrema or not.
  bool const cubic_has_local_extrema = d > 0.0;

  // Calculate the inflection point (where the second derivative is zero).
  double const Ix = -c2 / 3.0;
  // Magic (took me several days to find this).
  double const M = 27.0 * c0 + (2.0 * d - 3.0 * c1) * c2;

  double C0, C1;
  double scale;

  if (cubic_has_local_extrema)
  {
    // Transform the cubic into
    //
    //   Q(u) = C0 - 3u + u³
    //
    // After transforming the cubic we have the following four possibilities:
    //
    //          starting point
    //                ↓
    //  --------------+--O--> u-axis (for example)
    //      .-.         /
    //     /   \       /
    //    /     \←----/----- Inflection point
    //   /       \   /
    //            `-´←------ The extreme with the largest absolute y value.
    //       ↑  ↑  ↑  ↑
    //      -1  0  1  2
    //
    // The starting point is set on the positive slope on the side of the local extreme that is the furthest away from the x-axis.
    //
    double const sqrt_d = std::sqrt(d);
    C0 = M / (d * sqrt_d);
    C1 = -3.0;

    // The applied transform means that any root found must be scaled back by multiplying with
    scale = sqrt_d / 3.0;
    // and then adding Ix back.
  }
  else
  {
    // Transform the cubic into
    //
    //   Q(u) = C0 + 3u + u³
    //
    //          starting point
    //                ↓
    //                /
    //  -------------O+-> u-axis (for example)
    //              /
    //             /
    //            /
    //         .⋅´←--------- Inflection point
    //        /
    //       /
    //      /
    //       ↑  ↑  ↑  ↑
    //      -1  0  1  2

    double const sqrt_md = std::sqrt(-d);
    C0 = M / (-d * sqrt_md);
    C1 = 3.0;

    // The applied transform means that any root found must be scaled back by multiplying with
    scale = sqrt_md / 3.0;
    // and then adding Ix back.
  }

  // Determine if the inflection point is above or below the x-axis.
  bool const inflection_point_y_larger_than_zero = C0 > 0.0;

  // Special case for if zero is a root (i.e. Q(u) = u⋅(u² ± 3)).
  double u = 0.0;

  if (AI_LIKELY(C0 != 0.0))       // Is 0 not a root of the cubic?
  {
    // If cubic_has_local_extrema and C0 is less than 0.90375741845959156233304814223072905692
    // then the following fifth degree polynomial approximates S(C0) with a maximum absolute
    // error of 0.00012624.
    //
    // S_approximation(C0) =
    //   -0.0001262402246896549171879        +
    //    0.0605169674714465246944928 * C0   +
    //    1.3637323672416603782832849 * C0^2 +
    //   -1.1689405475539427807915498 * C0^3 +
    //    0.2645861875793937051035535 * C0^4 +
    //    0.0295306624154090842689822 * C0^5
    //
    // or, this third degree polynomial results in a maximum absolute error of 0.001783.
    //
    //   -0.0017830414728385209713872        +
    //    0.1217675925997833002784848 * C0   +
    //    1.0356078024797888438012922 * C0^2 +
    //   -0.6176223511426321157909154 * C0^3
    //
    // where S(C0) is defined as that (S(C0) - 1) sqrt(3) - S(C0) (cbrt(C0) + 1 / cbrt(C0)) is the exact root.

    // Avoid the local extrema and the inflection point because the derivative might be zero there too.
    double cbrtC0 = std::cbrt(C0);
    u = -cbrtC0 - 1.0 / cbrtC0;
    Dout(dc::notice, "Initial guess: " << u);

    // Since the initial guess is already very accurate, it is more than sufficient to
    // determine the resolution of a double around the value of the root (nextafter(u) - u).
    // If we add less than 1.5 times that to u (but more than 0.5 times that) then u will
    // be incremented with more than this resolution delta, which is when another iteration
    // can still improve the result.
    double const epsilon = 1.5 * (std::nextafter(u, std::numeric_limits<double>::infinity()) - u);
    ASSERT(epsilon > 0.0);
    int limit = 10000;
    double prev_u;
    double step = std::numeric_limits<double>::infinity();
    double prev_step;
    do
    {
      prev_u = u;
      prev_step = step;
      // Calculate Q(u) = C0 + C1 * u + u^3.
      double Q_u = C0 + u * (utils::square(u) + C1);
      // Calculate Q''(u) = 6 * u;
      double half_Qpp_u = 3.0 * u;
      // Calculate Q'(u) = C1 + 3 * u^2.
      double Qp_u = half_Qpp_u * u + C1;
      // Apply Halley's method.
      step = -Q_u * Qp_u / (utils::square(Qp_u) - Q_u * half_Qpp_u);
      u += step;                                                                // uₙ₊₁ = uₙ - Q(u)Q'(u) / (Q'(u)² - ½Q(u)Q"(u))
#ifdef CWDEBUG
      Dout(dc::notice, "Halley: u = " << std::setprecision(18) << u << " (" <<
          std::nextafter(u, std::numeric_limits<double>::infinity()) << "); step = " << step << "; Δu = " << (u - prev_u));
      // Make sure that comparing with epsilon doesn't do worse than detecting that u only changed by a single resolution delta.
      double near = step > 0.0 ? std::nextafter(prev_u, std::numeric_limits<double>::infinity())
                               : std::nextafter(prev_u, -std::numeric_limits<double>::infinity());
      // If u only changed a single bit, then step should NOT be larger than epsilon!
      if (u == near)
        ASSERT(!(std::abs(step) > epsilon));
#endif
    }
    while (std::abs(step) > epsilon && --limit);
    ASSERT(limit > 0);
    iterations = 10000 - limit;
  }

  roots_out[0] = u * scale + Ix;
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
