#pragma once

#include "utils/square.h"
#include "utils/macros.h"
#include <array>
#include <cmath>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class QuadraticPolynomial
{
 private:
  std::array<double, 3> coefficients_{};

 public:
  // Create a zero polynomial.
  QuadraticPolynomial() { }
  // Create a polynomial  a + b x + c x^2.
  QuadraticPolynomial(double a, double b, double c) : coefficients_{{a, b, c}} { }

  // Evaluation.
  double operator()(double w) const
  {
    return coefficients_[0] + (coefficients_[1] + coefficients_[2] * w) * w;
  };

  // Evaluate derivative.
  double derivative(double w) const
  {
    return coefficients_[1] + 2.0 * coefficients_[2] * w;
  }

  // Access coefficients.
  double operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  double& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

  int get_roots(std::array<double, 2>& roots_out) const
  {
    if (coefficients_[2] == 0.0)
    {
      roots_out[0] = -coefficients_[0] / coefficients_[1];
      return std::isfinite(roots_out[0]) ? 1 : 0;
    }

    double const D = utils::square(coefficients_[1]) - 4.0 * coefficients_[2] * coefficients_[0];
    if (D < 0.0)
      return 0;
    // Use a sqrt with the same sign as coefficients_[1];
    double const signed_sqrt_D = std::copysign(std::sqrt(D), coefficients_[1]);

    // Calculate the root closest to zero.
    roots_out[0] = -2.0 * coefficients_[0] / (coefficients_[1] + signed_sqrt_D);

    if (AI_UNLIKELY(std::isnan(roots_out[0])))
    {
      // This means we must have divided by zero, which means that both, coefficients_[1] as well as sqrtD, must be zero.
      // The latter means that coefficients_[0] is zero (coefficients_[2] was already checked not to be zero).
      // Therefore we have: f(x) = c x^2 with one root at x=0.
      roots_out[0] = 0.0;
      return 1;
    }

    // Calculate the root further away from zero.
    roots_out[1] = -0.5 * (coefficients_[1] + signed_sqrt_D) / coefficients_[2];

    // The second one is larger in absolute value.
    ASSERT(std::abs(roots_out[1]) > std::abs(roots_out[0]));

    return 2;
  }

  // Returns the x coordinate of the vertex of the parabola.
  double vertex_x() const
  {
    // f(x) = a + bx + cx^2
    // f'(x) = b + 2c x
    // f'(x) = 0 --> x = -b / 2c
    return -0.5 * coefficients_[1] / coefficients_[2];
  }

  // Returns the y coordinate of the vertex of the parabola.
  double vertex_y() const
  {
    // f(vertex_x()) = a + b(-b / 2c) + c(-b / 2c)^2 = a - b^2 / 2c + b^2 / 4c = a - b^2 / 4c.
    return coefficients_[0] - 0.25 * utils::square(coefficients_[1]) / coefficients_[2];
  }

  QuadraticPolynomial height() const
  {
    // f(x) - vertex_y()
    QuadraticPolynomial result(*this);
    result[0] = 0.25 * utils::square(coefficients_[1]) / coefficients_[2];
    return result;
  }

  QuadraticPolynomial& operator-=(QuadraticPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] -= rhs.coefficients_[i];
    return *this;
  }

  QuadraticPolynomial& operator+=(QuadraticPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] += rhs.coefficients_[i];
    return *this;
  }

  QuadraticPolynomial operator+(QuadraticPolynomial const& rhs) const
  {
    QuadraticPolynomial result(*this);
    result += rhs;
    return result;
  }

  QuadraticPolynomial operator-(QuadraticPolynomial const& rhs) const
  {
    QuadraticPolynomial result(*this);
    result -= rhs;
    return result;
  }

  QuadraticPolynomial& operator*=(double factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] *= factor;
    return *this;
  }

  friend QuadraticPolynomial operator*(double factor, QuadraticPolynomial const& rhs)
  {
    QuadraticPolynomial result(rhs);
    result *= factor;
    return result;
  }

  // Returns true if the old parabola is close enough to this parabola at infinity.
  // toggles_out is filled with an even number of x-coordinates at which the
  // y-coordinate of the old parabola goes in or out the acceptable range.
  bool equal_intervals(QuadraticPolynomial const& old, std::array<double, 4>& toggles_out, int& count_out) const
  {
    // This object is the most recent parabola, the one that came after old.

    // Create a new parabola, diff = 10 * (new - old).
    QuadraticPolynomial diff{10.0 * (coefficients_[0] - old[0]), 10.0 * (coefficients_[1] - old[1]), 10.0 * (coefficients_[2] - old[2])};

    // If abs(diff[2]) < abs(new[2]) then the old one is close enough at infinity.
    bool close_at_inf = std::abs(diff[2]) < std::abs(coefficients_[2]);

    // Set v_y to the y-coordinate of the vertex of the new parabola.
    double v_y = vertex_y();

    std::array<std::array<double, 2>, 2> toggle;
    std::array<int, 2> count{};
    std::array<double, 2> c{coefficients_[2] - diff.coefficients_[2], coefficients_[2] + diff.coefficients_[2]};
    std::array<double, 2> abs_c{std::abs(c[0]), std::abs(c[1])};
    for (int i = 0; i < 2; ++i) // i=0: subtract diff, i=1: add diff.
    {
      double sign = i == 0 ? -1.0 : 1.0;
      if (abs_c[i] > 1e-30)
      {
        double b = coefficients_[1] + sign * diff.coefficients_[1];
        double D = utils::square(b) - 4.0 * (coefficients_[0] - v_y + sign * diff.coefficients_[0]) * c[i];
        if (D > 0.0)
        {
          double avg = -0.5 * b / c[i];
          double delta = 0.5 * std::sqrt(D) / abs_c[i];
          toggle[i][0] = avg - delta;
          toggle[i][1] = avg + delta;
          count[i] = 2;
        }
      }
    }
    count_out = count[0] + count[1];
    int i0 = 0;
    int i1 = 0;
    int i = 0;

    // Sort the result and write it to toggles_out.
    while (i < count_out)
    {
      if (i1 == count[1] || (i0 < count[0] && toggle[0][i0] < toggle[1][i1]))
        toggles_out[i] = toggle[0][i0++];
      else
        toggles_out[i] = toggle[1][i1++];
      ++i;
    }

    ASSERT(count_out == 0 ||
           (count_out == 2 && toggles_out[0] < toggles_out[1]) ||
           (count_out == 4 && toggles_out[0] < toggles_out[1] && toggles_out[1] < toggles_out[2] && toggles_out[2] < toggles_out[3]));

    return close_at_inf;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
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
};

} // namespace math
