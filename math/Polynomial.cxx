#include "sys.h"
#include "Polynomial.h"
#include <Eigen/Dense>
#include "utils/macros.h"
#include "utils/square.h"

namespace math {

int Polynomial::get_roots(std::array<double, 2>& roots_out) const
{
  // This can be at most a parabola.
  ASSERT(1 <= coefficients_.size() && coefficients_.size() <= 3);
  if (coefficients_.size() < 3 || coefficients_[2] == 0.0)
  {
    if (coefficients_.size() < 2)
      return 0;
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

int Polynomial::get_roots(std::array<std::complex<double>, 5>& roots_out) const
{
  int degree = coefficients_.size() - 1;

  // This function only works up till degree 5.
  ASSERT(degree <= 5);

  // Abbreviate the coefficients as a[0] etc.
  auto const& a = coefficients_;

  if (degree == 5 && a[5] != 0.0)
  {
    // Construct the 5x5 companion matrix.
    Eigen::Matrix<double, 5, 5> C;
    C <<  0, 0, 0, 0, -a[0] / a[5],
          1, 0, 0, 0, -a[1] / a[5],
          0, 1, 0, 0, -a[2] / a[5],
          0, 0, 1, 0, -a[3] / a[5],
          0, 0, 0, 1, -a[4] / a[5];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix<double, 5, 5>> solver(C);
    Eigen::Matrix<std::complex<double>, 5, 1> roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
        roots_out[i] = roots[i];
  }
  else if (degree >= 4 && a[4] != 0.0)
  {
    // Construct the 4x4 companion matrix.
    Eigen::Matrix4d C;
    C <<  0, 0, 0, -a[0] / a[4],
          1, 0, 0, -a[1] / a[4],
          0, 1, 0, -a[2] / a[4],
          0, 0, 1, -a[3] / a[4];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix4d> solver(C);
    Eigen::Vector4cd roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
      roots_out[i] = roots[i];
  }
  else if (degree >= 3 && a[3] != 0.0)
  {
    degree = 3;

    // Construct the 3x3 companion matrix.
    Eigen::Matrix3d C;
    C <<  0, 0, -a[0] / a[3],
          1, 0, -a[1] / a[3],
          0, 1, -a[2] / a[3];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix3d> solver(C);
    Eigen::Vector3cd roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
      roots_out[i] = roots[i];
  }
  else if (degree >= 2 && a[2] != 0.0)
  {
    degree = 2;

    // Construct the 2x2 companion matrix.
    Eigen::Matrix2d C;
    C <<  0, -a[0] / a[2],
          1, -a[1] / a[2];

    // Compute the eigenvalues of the companion matrix.
    Eigen::EigenSolver<Eigen::Matrix2d> solver(C);
    Eigen::Vector2cd roots = solver.eigenvalues();

    // Store the roots in the output array.
    for (int i = 0; i < roots.size(); ++i)
      roots_out[i] = roots[i];
  }
  else if (degree >= 1 && a[1] != 0.0)
  {
    degree = 1;

    roots_out[0] = -a[0] / a[1];
  }
  else
    degree = 0;

  return degree;
}

#ifdef CWDEBUG
void Polynomial::print_on(std::ostream& os) const
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
        os << ' ' << symbol_name_;
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
