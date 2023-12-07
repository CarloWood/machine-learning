#include "gnuplot_i.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <array>
#include <sstream>

//
//                      N points
// <------------------------------------------------------->
// + + + + + + + + + + + + + + + + +.+.+.+.+.+.+.+.+.+.+.+.+ ^
// + + + + + + + + + + + + + + +.+ + + + + + + + + + + + + + |
// + + + + + + + + + + + + + +.+ + + + + + + + + + + + + + + |
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + |
// + + + + + + + + + + + + +.+ + + + + + + + + + + + + + + + |M points
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + |
// + + + + + + + + + + + +.+ + + + + + + + + + + + + + + + + |
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + |
// + + + + + + + + + + +.+ + + + + + + + + + + + + + + + + + |
// + + + + + + + + + +.+ + + + + + + + + + + + + + + + + + + |
// +.+.+.+.+.+.+.+.+ + + + + + + + + + + + + + + + + + + + + v
// -------------------------------------------------------------> x
// ^         | |      <----|----> 2σ = 2sπ/√3              ^
// |0         dx           mu                              |total_width
// <----------------------------- Sσ ---------------------->
//
constexpr int M = 10;
constexpr long double S = 40;           // Integrate over +/- 20 standard deviations.
constexpr long double s = 2.0;
constexpr long double pidsqrt3 = 1.813799364234217850594078258L; // pi / sqrt(3)
constexpr long double sd = s * pidsqrt3;
constexpr long double total_width = S * sd;
constexpr long double dx = 0.05L * s; //total_width / (N - 1);
constexpr int N = 1 + total_width / dx;
constexpr long double mu = total_width * 0.5;
constexpr long double wa = 1 / s;
constexpr long double ba = -mu / s;     // b = -β/s, and β=μ in this case because γ=0.
constexpr long double pi2d6 = 1.6449340668482264364724151666460L; // pi^2 / 6
constexpr long double epsilon = 1e-14L; // Can be used as abs_relative_error.

// Copied from utils/almost_equal.h
//
// Return true when
//
//      |z1 - z2| <= (|z1 + z2| / 2) * abs_relative_error.
//
// In other words, the flipping point occurs when
//                                                         |    z1 - z2    |
//      abs_relative_error = |z1 - z2| / (|z1 + z2| / 2) = | ------------- |
//                                                         | (z1 + z2) / 2 |
//
template<class T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type
   almost_equal(T x, T y, T const abs_relative_error)
{
  return 2 * std::abs(x - y) <= abs_relative_error * std::abs(x + y);
}

long double sigmoid(long double x)
{
  return 1.0L / (1.0L + std::exp(-x));
}

// Calculate
//    ∞
//   ⌠
//   ⎮ σ(u) Log(σ(weight_ratio u + bias_error)) du
//   ⌡
// -∞
// Where weight_ratio = w₁/wₐ and bias_error = b - bₐ w₁/wₐ.
//
long double integral_sigma_log_sigma(long double weight_ratio, long double bias_error)
{
  long double eps = epsilon / N;
  // u  = x wa + ba
  // du = wa dx
  constexpr long double du = wa * dx;
  long double sum = 0;
  for (int i = 0;; ++i)
  {
    long double u = i * du;
    long double ds = sigmoid(u) * std::log(sigmoid(u * weight_ratio + bias_error));
    sum += ds;
    if (i >= N && std::abs(ds) < eps)
      break;
  }
  for (int i = -1;; --i)
  {
    long double u = i * du;
    long double ds = sigmoid(u) * std::log(sigmoid(u * weight_ratio + bias_error));
    sum += ds;
    if (i <= -N && std::abs(ds) < eps)
      break;
  }
  return sum * du;
}

// The same, for bias_error = 0.
//
// Integral from -inf to +inf of [sigmoid(u) log(sigmoid(n u)) du] = -PI^2 a(n) / 24n.
// Also works for non-integer n>0, when a(n)=2n^2+2. See https://oeis.org/draft/A005893
long double integral_sigma_log_sigma(long double weight_ratio)
{
  long double A005893 = 2 * (weight_ratio * weight_ratio + 1);
  long double denominator = -4.0L * weight_ratio;
  return A005893 * pi2d6 / denominator;
}

long double g_integral(long double weight_ratio, long double bias_error)
{
  long double g_weight_ratio_0 = integral_sigma_log_sigma(weight_ratio);
  return g_weight_ratio_0 - 0.25L / weight_ratio * bias_error * bias_error;
}

int main()
{
  std::cout.precision(12);

  std::cout << "N = " << N << std::endl;
  std::cout << "M = " << M << std::endl;
  std::cout << "s = " << s << " (sd = " << sd << ")" << std::endl;
  std::cout << "total_width = " << total_width << std::endl;
  std::cout << "mu = " << mu << std::endl;
  std::cout << "dx = " << dx << std::endl;

#if 0   // Verify that we can numerically calculate the integral with enough precision.
  long double sum = 0;
  for (int i = 0; i < N; ++i)
  {
    long double x = i * dx;
    long double ds = sigmoid(x * wa + ba) * std::log(sigmoid(x * wa + ba));
    sum += ds * dx;
  }
  std::cout << "sum / s = " << (sum / s) << std::endl;
  std::cout << "sum  + s * PI^2 / 6 = " << (sum + s * pi2d6) << std::endl;
  if (std::abs(sum + s * pi2d6) < 1e-12L)
    std::cout << "Precision of 12 digits reached." << std::endl;
  assert(almost_equal(sum, -s * pi2d6, epsilon));
#endif

#if 0   // Calculate Log(likeliness_opt) - see README.binary_classification
  long double log_likeliness_opt = 0.0L;
  for (int i = 0; i < N; ++i)
  {
    // x = i Δx
    long double x = i * dx; // x runs from 0 till total_width.
//    std::cout << "x = " << x << std::endl;

    // P_greenᵢ = σ(x w₁ + b), here with w₁ = wₐ and b = bₐ.
    long double P_green_i = sigmoid(x * wa + ba);
//    std::cout << "P_green_i = " << P_green_i << std::endl;

    // N_greenᵢ = M σ(x wₐ + bₐ)
    long double N_green_i = M * P_green_i;
    // P_redᵢ = 1 - P_greenᵢ
    long double P_red_i = 1.0L - P_green_i;
    // N_redᵢ = M - N_greenᵢ
    long double N_red_i = M - N_green_i;
    // likelinessᵢ = P_greenᵢ^N_greenᵢ * P_redᵢ^N_redᵢ
    long double log_likeliness_opt_i = std::log(P_green_i) * N_green_i + std::log(P_red_i) * N_red_i;
//    std::cout << "log_likeliness_opt_i(i=" << i << ") = " << log_likeliness_opt_i << std::endl;

    // Log(likeliness) = 𝚺 Log(likelinessᵢ)
    //                   i
    log_likeliness_opt += log_likeliness_opt_i;
  }
  //                                   ∞
  //                                  ⌠                         π² M
  // Log(likeliness_opt) = 2M/(Δx wₐ) ⎮ σ(u) Log(σ(u)) du = - ⎼⎼⎼⎼⎼⎼⎼⎼
  //                                  ⌡                        3 Δx wₐ
  //                                -∞

  // This prints -π²/3
  std::cout << "log(likeliness_opt) / (M / Δx wₐ) = " << (log_likeliness_opt / (M / (dx * wa))) << std::endl;
  assert(almost_equal(log_likeliness_opt, -2.0L * pi2d6 * M / (dx * wa), epsilon));
#endif

#if 0   // Verify integral_sigma_log_sigma gives the same result as the numerically calculated integral.
  for (long double we = 0.1L; we < 10.0L; we *= 1.1L)
  {
    // This should print the same value twice.
    double long int0 = integral_sigma_log_sigma(we, 0.0L);
    std::cout << integral_sigma_log_sigma(we, 0.0L) << std::endl;
    double long exact = integral_sigma_log_sigma(we);
    std::cout << exact << std::endl;
    assert(almost_equal(int0, exact, epsilon));
  }
#endif

#if 0   // Comparing g(1, be) with quadratic guess.
  std::vector<double> x;
  std::vector<double> y1;
  for (long double be = -100.0L; be < 100.01L; be += 0.1L)
  {
    long double first_int = integral_sigma_log_sigma(1.0L, be);
    long double second_int = integral_sigma_log_sigma(1.0L, -be);
    long double g1be = (first_int + second_int) / 2;

    std::cout << (-g1be / pi2d6) << std::endl;
    x.push_back(be);
    y1.push_back(be * be / 4 + g1be + pi2d6);
    // The guess is correct:
    assert(almost_equal(g1be, -pi2d6 - be * be / 4, epsilon));
  }
  gnuplot_ctrl* h1 = gnuplot_init();
  gnuplot_setstyle(h1, "lines");
  gnuplot_plot_coordinates(h1, x.data(), y1.data(), x.size(), "b factor");

  // Close plot window.
  std::cin.get();
  gnuplot_close(h1);
#endif

#if 1   // Comparing g(2, bias_error) with quadratic guess.
  std::vector<double> x;
  std::array<std::vector<double>, 4> ys;
  for (int wo = 0; wo < ys.size(); ++wo)
  {
    long double weight_ratio = 0.0625L * (1 << (2 * wo));       // Test with w₁/wₐ = 1/16, 1/4, 1 and 4.
    for (long double bias_error = -4.0L; bias_error < 4.01L; bias_error += 0.01L)
    {
      long double first_int = integral_sigma_log_sigma(weight_ratio, bias_error);
      long double second_int = integral_sigma_log_sigma(weight_ratio, -bias_error);
      long double g_weight_ratio_bias_error = (first_int + second_int) / 2;
      long double g = g_integral(weight_ratio, bias_error);

      std::cout << (-g_weight_ratio_bias_error / pi2d6) << std::endl;
      if (wo == 0)
        x.push_back(bias_error);
      ys[wo].push_back(g_weight_ratio_bias_error - g);
      // The guess is correct:
      assert(almost_equal(g_weight_ratio_bias_error, g, epsilon));
    }
  }
  gnuplot_ctrl* h1 = gnuplot_init();
  gnuplot_setstyle(h1, "linepoints");
  for (int wo = 0; wo < ys.size(); ++wo)
  {
    std::ostringstream title;
    long double weight_ratio = 0.0625L * (1 << (2 * wo));
    title << "wr=" << weight_ratio << std::endl;
    gnuplot_plot_coordinates(h1, x.data(), ys[wo].data(), x.size(), title.str().c_str());
  }

  // Close plot window.
  std::cin.get();
  gnuplot_close(h1);
#endif
}
