#pragma once

#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "utils/square.h"
#include "utils/macros.h"
#include <array>
#include <cmath>
#include <cstring>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class CubicPolynomial
{
 protected:
  std::array<double, 4> coefficients_{};

 public:
  // Create a zero polynomial.
  CubicPolynomial() { }
  // Create a polynomial  a + b x + c x^2 + d x^3.
  CubicPolynomial(double a, double b, double c, double d) : coefficients_{{a, b, c, d}} { }

  // Construct a cubic polynomial from a Polynomial.
  CubicPolynomial(Polynomial const& polynomial)
  {
    std::vector<double> const& coefficients = polynomial.coefficients();
    // polynomial must have degree three or less.
    ASSERT(coefficients.size() == 4);
    std::memcpy(coefficients_.data(), coefficients.data(), coefficients_.size() * sizeof(double));
  }

  void initialize(double x0, double y0, double dxdy0, double x1, double y1, double dxdy1);
  int get_extrema(std::array<double, 2>& extrema_out, bool left_most_first = true) const;

  // Return the derivative of this cubic.
  QuadraticPolynomial derivative() const
  {
    return {coefficients_[1], 2.0 * coefficients_[2], 3.0 * coefficients_[3]};
  }

  // Evaluation.
  double operator()(double x) const
  {
    return coefficients_[0] + (coefficients_[1] + (coefficients_[2] + coefficients_[3] * x) * x) * x;
  };

  // Evaluate derivative.
  double derivative(double x) const
  {
    return coefficients_[1] + (2.0 * coefficients_[2] + 3.0 * coefficients_[3] * x) * x;
  }

  // Evaluate second derivative.
  double second_derivative(double x) const
  {
    return 2.0 * coefficients_[2] + 6.0 * coefficients_[3] * x;
  }

  double inflection_point() const
  {
    if (AI_UNLIKELY(coefficients_[3] == 0.0))
      return std::numeric_limits<double>::infinity();           // Anything really large would do (or towards negative infinity too).
    return -coefficients_[2] / (3.0 * coefficients_[3]);
  }

  // Access coefficients.
  double operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  double& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

  CubicPolynomial& operator-=(CubicPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] -= rhs.coefficients_[i];
    return *this;
  }

  CubicPolynomial& operator+=(CubicPolynomial const& rhs)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] += rhs.coefficients_[i];
    return *this;
  }

  CubicPolynomial operator+(CubicPolynomial const& rhs) const
  {
    CubicPolynomial result(*this);
    result += rhs;
    return result;
  }

  CubicPolynomial operator-(CubicPolynomial const& rhs) const
  {
    CubicPolynomial result(*this);
    result -= rhs;
    return result;
  }

  CubicPolynomial& operator*=(double factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] *= factor;
    return *this;
  }

  friend CubicPolynomial operator*(double factor, CubicPolynomial const& rhs)
  {
    CubicPolynomial result(rhs);
    result *= factor;
    return result;
  }

  // Return the division of this Polynomial by the factor (x - r).
  QuadraticPolynomial long_division(double r, double& remainder) const
  {
    QuadraticPolynomial result;
    result[2] = coefficients_[3];
    result[1] = coefficients_[2] + r * result[2];
    result[0] = coefficients_[1] + r * result[1];
    remainder = coefficients_[0] + r * result[0];
    return result;
  }

  int get_roots(std::array<double, 3>& roots_out, int& iterations) const;
  int get_roots(std::array<double, 3>& roots_out) const
  {
    int iterations;
    return get_roots(roots_out, iterations);
  }

#if CW_DEBUG
  // Return true if one was assigned from the other.
  friend bool operator==(CubicPolynomial const& lhs, CubicPolynomial const& rhs)
  {
    return lhs.coefficients_ == rhs.coefficients_;
  }

  void print_on(std::ostream& os) const;
#endif
};

} // namespace math
