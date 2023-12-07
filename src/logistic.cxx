#include "gnuplot_i.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <random>
#include <thread>
#include <chrono>
#include <fstream>

double const one_over_sqrt_of_two_pi = 1.0 / std::sqrt(2.0 * M_PI);

double normal_distribution_0_1(double x)
{
  double minus_half_x2 = -0.5 * x * x;
  return one_over_sqrt_of_two_pi * std::exp(minus_half_x2);
}

double logistic_distribution_0(double x, double s)
{
  double emx = std::exp(-x/s);
  double opemx = 1.0 + emx;
  return emx / (s * opemx * opemx);
}

double sigmoid(double x)
{
  double emx = std::exp(-x);
  double omemx = 1.0 + emx;
  return 1.0 / omemx;
}

constexpr double inf = 100.0;
constexpr double dx = 0.001;

double integral(double from, double to, std::function<double(double)> f)
{
  if (from < -100.0)
    from = -inf;
  if (to > 100.0)
    to = inf;
  double sum = 0;
  for (double x = from; x < to; x += dx)
    sum += f(x) * dx;
  return sum;
}

void draw_points(gnuplot_ctrl* h1, std::array<std::vector<double>, 2> const& xp, std::array<std::vector<double>, 2> const& yp)
{
  gnuplot_setstyle(h1, "points");
  gnuplot_append_style(h1, " pointtype 7 pointsize 2 linewidth 2 linecolor 'red'");
  gnuplot_plot_coordinates(h1, xp[0].data(), yp[0].data(), xp[0].size(), "inputs0");

  gnuplot_setstyle(h1, "points");
  gnuplot_append_style(h1, " pointtype 4 pointsize 2 linewidth 2 linecolor 'green'");
  gnuplot_plot_coordinates(h1, xp[1].data(), yp[1].data(), xp[1].size(), "inputs1");

  gnuplot_setstyle(h1, "points");
  gnuplot_cmd(h1, "unset key");
//  gnuplot_cmd(h1, "set xrange [0:10]");
  gnuplot_cmd(h1, "set yrange [0:1]");
}

void draw_function(gnuplot_ctrl* h1, std::string const& function_name, std::string const& equation)
{
  gnuplot_cmd(h1, equation.c_str());
  gnuplot_setstyle(h1, "lines");
  gnuplot_plot_equation(h1, function_name.c_str(), equation.c_str());
}

double log_likeliness(double mu, double s, std::vector<std::pair<double, bool>> const& points)
{
  double result = 0.0L;
  for (auto const& point : points)
  {
    double x = point.first;
    bool is_green = point.second;
    double vx = is_green ? x - mu : mu - x;
    result += std::log(sigmoid(vx / s));
  }
  return result;
}

int main()
{
  std::cout.precision(9);

  // GNU plot handle.
  gnuplot_ctrl* h1 = gnuplot_init();

#if 0
  // Normal distribution function.
  {
    double const sd = 1;
    double variance = 0;
    for (double x = -100.0; x < 100.0; x += dx)
    {
      double nd = normal_distribution_0_1(x);
      variance += x * x * nd * dx;
    }
    double standard_deviation = std::sqrt(variance);
    std::cout << "Standard deviation = " << standard_deviation << std::endl;
    std::cout << "sd = " << sd << std::endl;
    double cdf = 0;
    for (double x = -100.0; x < sd; x += dx)
    {
      double nd = normal_distribution_0_1(x);
      cdf += nd * dx;
    }
    std::cout << "Integral -inf till sd of pdf(x) dx = " << cdf << std::endl;
    std::cout << "0.5 + 34.1 / 100 = 0.841" << std::endl;
  }
#endif

#if 1
  // Logistic distribution function.
  {
    for (double s = 0.25; s < 5.0; s *= 2)
    {
      double sd = s * M_PI / std::sqrt(3);
      double Psd1 = integral(-inf, sd, [s](double x){ return logistic_distribution_0(x, s); });
      double Psd2 = integral(-inf, sd / s, [s](double x){ return logistic_distribution_0(x, 1.0); });
      std::cout << Psd1 << " : " << Psd2 << std::endl;
    }
    std::cout << sigmoid(M_PI / std::sqrt(3)) << std::endl;

    for (double s = 0.25; s < 5.0; s *= 2)
    {
      double variance = integral(-inf, inf, [s](double x){ return x * x * logistic_distribution_0(x, s); });
      std::cout << "sqrt(variance) / s = " << (std::sqrt(variance) / s) << std::endl;
    }

    for (double s = 0.25; s < 5.0; s *= 2)
    {
      double average_distance = integral(-inf, inf, [s](double x){ return std::abs(x) * logistic_distribution_0(x, s); });
      std::cout << "average_distance / s = " << (average_distance / s) << std::endl;
    }
  }
#endif

  // Generate N data points between 0 and 10.
  constexpr int N = 1000;
  std::vector<std::pair<double, bool>> points;

  std::random_device rd;
  std::default_random_engine generator(rd());

  // Put the mean mu between 3 and 7.
  std::uniform_real_distribution<double> mu_dist(3.0, 7.0);
  double mu = mu_dist(generator);
  mu = 0;
  std::cout << "mu = " << mu << std::endl;

  // Generate the scaling factor between 0.5 and 2.
  std::uniform_real_distribution<double> s_dist(-1.0, 2.0);
  double s = std::exp(std::log(2) * s_dist(generator));
  s = 1;
  double sd = s * M_PI / std::sqrt(3);
  std::cout << "Standard deviation = " << sd << std::endl;

  // All inputs x are in the range mu +/- 5 * sigma.
  std::uniform_real_distribution<double> x_dist(mu - 5 * sd, mu + 5 * sd);
  // Sigmoid returns value between 0 and 1. We need that too.
  std::uniform_real_distribution<double> color_dist(0.0, 1.0);

  // X and Y coordinates to plot red (index = 0) and green (index = 1) points.
  std::array<std::vector<double>, 2> xp;
  std::array<std::vector<double>, 2> yp;

  for (int i = 0; i < N; ++i)
  {
    // Pick an x.
    double x = x_dist(generator);
    // Calculate the probability that this point is green.
    double Pgreen = sigmoid((x - mu) / s);
    // Determine if the point is red or green.
    double color = color_dist(generator);
    bool is_green = color < Pgreen;

    // Store all generated points.
    points.emplace_back(x, is_green);

    xp[is_green].push_back(x);
    yp[is_green].push_back(color);
  }
#if 1
  draw_points(h1, xp, yp);
  std::ostringstream equation;
  equation << "f(x) = 1 / (1 + exp(-(x - " << mu << ")/" << s << "))";
  draw_function(h1, "f(x)", equation.str());
  //std::this_thread::sleep_for(std::chrono::milliseconds(100));
  //gnuplot_resetplot(h1);
#endif

  long double expected_best_log_likeliness = log_likeliness(mu, s, points);
  std::cout << "ln(best_likeliness) (s = " << s << ") = " << expected_best_log_likeliness << std::endl;

  gnuplot_ctrl* h2 = gnuplot_init();
  std::vector<double> mu_grid;
  std::vector<double> s_grid;
  std::vector<double> likeliness;
  double delta_s = 0.2 / std::sqrt(N);
  double delta_mu = 0.2 / std::sqrt(N);
  // In delta_s' offsets relative to s.
  int s_start = -50;
  int s_end = 49;
  // In delta_mu's offsets relative to mu.
  int mu_start = -50;
  int mu_end = 49;
  int const rows = mu_end - mu_start + 1;
  int const cols = s_end - s_start + 1;
  double test_mu = mu + mu_start * delta_mu;
  for (int row = 0; row < rows; test_mu += delta_mu, ++row)
  {
    double test_s = s + s_start * delta_s;
    for (int col = 0; col < cols; test_s += delta_s, ++col)
    {
      double loss = expected_best_log_likeliness - log_likeliness(test_mu, test_s, points);
//      std::cout << "ln(likeliness / best_likeliness) (mu = " << test_mu << "; s = " << test_s << ") = \t" << loss << std::endl;
      mu_grid.push_back(test_mu);
      s_grid.push_back(test_s);
      if (loss > 4)
        loss = 4;
      likeliness.push_back(loss);
    }
  }
  // Find the minimum.
  int mi = 0;
  double ml = likeliness[mi];
  for (int i = 0; i < (int)likeliness.size(); ++i)
    if (likeliness[i] < ml)
    {
      ml = likeliness[i];
      mi = i;
    }
  double best_mu = mu_grid[mi];
  double best_s = s_grid[mi];
  int mi_row = mi / cols;
  int mi_col = mi % cols;
  // Prepare data for likeliness(s) through best mu value.
  std::vector<double> s_best_mu;
  std::vector<double> likeliness_best_mu;
  double last_loss = 0;
  double max_slope = 0;
  for (double test_s = 0.01; test_s <= s + 5 * sd; test_s += 0.1 * delta_s)
  {
    s_best_mu.push_back(test_s);
    double loss = -log_likeliness(best_mu, test_s, points) * test_s;
    if (test_s == 0.01)
      std::cout << "loss(s) = N * " << (loss / N) << " / s, around s=0.\n";
//    if (test_s > 1.0)
    {
      if (last_loss != 0)
      {
        double slope = (loss - last_loss) / (0.1 * delta_s);
//        likeliness_best_mu.push_back(slope);
        if (slope > max_slope)
          max_slope = slope;
      }
//      else
//        likeliness_best_mu.push_back(0);
      last_loss = loss;
    }
    likeliness_best_mu.push_back(loss / N);
  }
  std::cout << "slope = N * " << (max_slope / N) << std::endl;
  gnuplot_ctrl* h3 = gnuplot_init();
  gnuplot_setstyle(h3, "lines");
  gnuplot_plot_coordinates(h3, s_best_mu.data(), likeliness_best_mu.data(), s_best_mu.size(), "loss(s|best mu)s/N");

#if 0
  bool plot_3d = rows > 1 && cols > 1;
  if (plot_3d)
  {
    //gnuplot_cmd(h2, "set view equal xyz");
    gnuplot_cmd(h2, "set view 60, 30");

    bool plot_grid = false;
    if (plot_grid)
    {
      gnuplot_cmd(h2, "set hidden3d");
      gnuplot_cmd(h2, "set pm3d");
      gnuplot_cmd(h2, "set cbrange [-1:1]");
      gnuplot_setstyle(h2, "lines");
      gnuplot_splot_grid(h2, likeliness.data(), rows, cols, "likeliness");
    }
    else
    {
      gnuplot_setstyle(h2, "points");
      gnuplot_splot(h2, mu_grid.data(), s_grid.data(), likeliness.data(), rows * cols, "likeliness");
    }
  }
  else
  {
    gnuplot_cmd(h2, "set xrange [0.9:1.1]");
    gnuplot_setstyle(h2, "lines");
    gnuplot_plot_coordinates(h2, (rows > 1) ? mu_grid.data() : s_grid.data(), likeliness.data(), rows * cols, "likeliness");
  }
#endif

  // Close plot window.
  std::cin.get();
  gnuplot_close(h3);
  gnuplot_close(h2);
  gnuplot_close(h1);
}
