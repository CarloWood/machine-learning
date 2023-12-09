#include "intersection_points.h"
#include "gnuplot_i.h"
#include "tfdebug/DebugTensor.h"
#include <cmath>
#include <random>
#include <utility>
#include <iostream>
#include <sstream>
#include <cassert>

// This program tests the loss function as found in README.binary_classification.

using float_type = long double;
constexpr float_type Pi = 3.1415926535897932384626433832795029L;

using RandomGenerator = std::mt19937;

class Color
{
 private:
  bool is_green_;

 public:
  constexpr Color(bool is_green) : is_green_(is_green) { }

  bool is_green() const
  {
    return is_green_;
  }
};

constexpr Color red(false);
constexpr Color green(true);

// The standard deviation of a logistic distribution with scale factor s=1.
float_type const sd1 = Pi / std::sqrt(3.0L);

// Define the logistic sigmoid function.
float_type sigmoid(float_type u)
{
  return 1 / (1 + std::exp(-u));
}

// The line that devides red and green.
//
// The line exists of points (x₀, x₁) on a straight line.
//
// We define it here by means of a vector W = [w₀, w₁, w₂] such that
// W·X = 0 (the dot-product), where X = [x₀, x₁, 1].
//
// W·X = w₀x₀ + w₁x₁ + w₂ = 0 --> x₁ = (-w₀/w₁)x₀ - w₂/w₁
//
// If we define β = -w₂/w₁ and γ = -w₀/w₁ then the line
// could be written as x₁ = γ x₀ + β. However, then the line
// wouldn't be defined if w₁=0, and calculations would get
// unstable when w₁ is close to zero. Using W·X = 0 as the
// definition of the line does not lead to such problems.
//
// Notice that scaling W doesn't change the line.
//
// If we want to generate a random line with an equal chance
// for any angle θ that it makes with the positive-x axis
// then first we should generate a random θ in the range [0, 2π>
// and then set at first w₀=-sin(θ) and w₁=cos(θ).
// The unit normal to the line is then N = (w₀, w₁).
//
// Let dN be on the line, where |d| is the distance from
// the origin to the line. Then w₀(hw₀) + w₁(hw₁) + w₂ = 0 -->
// h(w₀²+w₁²) + w₂ = 0 --> h = -w₂.
//
// As this is with w₀²+w₁² = 1, we have as scaling factor
// s = 1/sqrt(w₀² + w₁²) = 1. And the standard deviation
// that belongs to the blurring is sd=π/√3.
//
class Line
{
 private:
  Tensor<2> W_;    // A row vector of three elements.

 public:
  Line(RandomGenerator& generator) : W_({1, 3})
  {
    // The orientation of the line has a uniformly distributed angle.
    std::uniform_real_distribution<float_type> theta_dist(-0.5L * Pi, 1.5L * Pi);
    float_type theta = theta_dist(generator);                     // Random angle between -π/2 and 3π/2.

    // Generate a scale factor between 0.25 and 4.
    std::uniform_real_distribution<float_type> scaling_factor_dist(-2, 2);
    float_type s = std::pow(2, scaling_factor_dist(generator));  // Random scaling factor between 2^-2 and 2^2.

    // Generate lines that have a random distance to the origin in the range [-2sd, 2sd].
    std::uniform_real_distribution<float_type> distance_dist(-2 * sd1 * s, 2 * sd1 * s);
    float_type distance = distance_dist(generator);

    W_(0, 0) = -sin(theta) / s;
    W_(0, 1) = cos(theta) / s;
    W_(0, 2) = -distance / s;           // Shift the line along the normal plus or minus twice the standard deviation.

    // N = (w0 * s, w1 * s)
    // kN = (k * w0 * s, k * w1 * s)
    // W kN = k * w0^2 * s + k * w1^2 * s = (w0^2 + w1^2) * k * s = k / s
    // W kN + w2 = 0 --> k / s + w2 = 0 --> k = -w2 * s = distance.

    std::cout << "s = " << s << "; distance = " << distance << std::endl;
  }

  float_type scale_factor() const
  {
    return 1 / std::sqrt(W_(0, 0) * W_(0, 0) + W_(0, 1) * W_(0, 1));
  }

  // Return W as a row vector.
  Tensor<2> const& W() const
  {
    return W_;
  }

  // Return a normal vector with length 1/s.
  Vector N() const
  {
    Vector N({3});
    N(0) = W_(0, 0);
    N(1) = W_(0, 1);
    N(2) = 0;
    return N;
  }

  Vector xrange() const
  {
    Vector result({2});
    float_type s = scale_factor();
    result(0) = -3 * sd1 * s;
    result(1) = 3 * sd1 * s;
    return result;
  }

  Vector yrange() const
  {
    Vector result({2});
    float_type s = scale_factor();
    result(0) = -3 * sd1 * s;
    result(1) = 3 * sd1 * s;
    return result;
  }

  bool draw(gnuplot_ctrl* h1, double x_min, double y_min, double x_max, double y_max) const;

  void print_on(std::ostream& os) const
  {
    float_type s = scale_factor();
    os << "{W_:" << W_ << "; s = " << s << "; sd = " << (sd1 * s) << '}';
  }

  friend std::ostream& operator<<(std::ostream& os, Line const& line)
  {
    line.print_on(os);
    return os;
  }
};

bool Line::draw(gnuplot_ctrl* h1, double x_min, double y_min, double x_max, double y_max) const
{
  float_type s = scale_factor();
  for (int sd = -1; sd <= 1; ++sd)
  {
    intersections::HyperPlane<float_type, 2> line({W_(0, 0), W_(0, 1)}, W_(0, 2) + sd * sd1);
    intersections::HyperBlock<float_type, 2> rectangle({x_min, y_min}, {x_max, y_max});

    auto intersections = rectangle.intersection_points(line);
    if (intersections.size() < 2)
      return false;
    if (intersections.size() != 2)
    {
      for (int i = 0; i < (int)intersections.size(); ++i)
        std::cout << intersections[i][0] << ", " << intersections[i][1] << std::endl;
    }

    {
      std::ostringstream oss;
      oss << "set arrow from " << intersections[0][0] << ", " << intersections[0][1] <<
        " to " << intersections[1][0] << ", " << intersections[1][1] << " nohead";
      gnuplot_cmd(h1, oss.str().c_str());
    }
  }
  gnuplot_cmd(h1, "set style line 1 lt 1 lc rgb \"purple\" lw 1");
  gnuplot_cmd(h1, "set grid back ls 1");
  gnuplot_cmd(h1, "set xtics 1");
  gnuplot_cmd(h1, "set ytics 1");
  {
    std::ostringstream oss;
    double distance = std::abs(W_(0, 2) * s);
    oss << "set object 1 circle at 0,0 size " << distance << " fc rgb \"purple\" lw 1";
    gnuplot_cmd(h1, oss.str().c_str());
  }
  gnuplot_cmd(h1, "set size ratio -1");
  gnuplot_cmd(h1, "replot NaN");

  return true;
}

std::vector<Vector> generate_points(RandomGenerator& generator, Line const& line, int number_of_points)
{
  // Push back number_of_points times (0, 0, 1) onto the points vector.
  Vector zero({3});
  zero(2) = 1;
  std::vector<Vector> points(number_of_points, zero);

  // Put all points evenly distributed in a square.
  Vector xrange = line.xrange();
  Vector yrange = line.yrange();
  std::uniform_real_distribution<float_type> x_position_dist(xrange(0), xrange(1));
  std::uniform_real_distribution<float_type> y_position_dist(yrange(0), yrange(1));

  for (int i = 0; i < number_of_points; ++i)
  {
    points[i](0) = x_position_dist(generator);
    points[i](1) = y_position_dist(generator);
  }

  return points;
}

std::vector<Color> generate_colors(RandomGenerator& generator, Line const& line, std::vector<Vector> const& points)
{
  std::vector<Color> targets;
  float_type s = line.scale_factor();
  float_type s3_inverse = 1 / (s * s * s);

  // Random number between 0 and 1.
  std::uniform_real_distribution<float_type> color_dist(0, 1);

  // Let W·X = 0 be the line equation where W = (w₀, w₁, w₂).
  // Let N = (w₀, w₁, 0) be the normal vector with length 1/s.
  // Let point be P = (p₀, p₁, 1).
  // Let K = P + kN be on the line, i.e. W·K = 0, then k is the "distance"
  // from point P to the line in units of N (the actual distance is |kN| = |k|/s).
  // Thus W·(P + kN) = 0 --> W·P + kW·N = 0 --> k = -W·P/W·N.
  // Note that W·N = w₀² + w₁² + w₂·0 = 1/s².
  // Hence, the real signed distance is k/s = -W·P s.
  // Finally, we have to devide the real distance by s before passing it to the sigmoid function.
  //
  // Calculate the scalar WN:
  Tensor W = line.W();                          // Row vector with shape {1,3}.
//  auto yr = line.yrange();
  for (Vector const& P : points)
  {
    float_type WP = contract(W, P, 1, 0)(0);    // W is a row vector, we have to run over its columns.
    float_type threshold = sigmoid(-WP);
    targets.push_back((color_dist(generator) < threshold) ? green : red);
//    targets.push_back(((P(1) - yr(0))/(yr(1) - yr(0)) < threshold) ? green : red);
  }

  return targets;
}

void plot(Line const& line, std::vector<Vector> const& points, std::vector<Color> const& targets)
{
  assert(points.size() == targets.size());
  std::array<std::vector<double>, 2> xp;
  std::array<std::vector<double>, 2> yp;
  for (int i = 0; i < points.size(); ++i)
  {
    Vector point = points[i];
    int color = targets[i].is_green() ? 1 : 0;
    xp[color].push_back(point(0));
    yp[color].push_back(point(1));
  }

  // GNU plot handle.
  gnuplot_ctrl* h1 = gnuplot_init();
  gnuplot_cmd(h1, "unset key");

  // Draw the red points.
  gnuplot_setstyle(h1, "points");
  gnuplot_append_style(h1, " pointtype 7 pointsize 2 linewidth 2 linecolor 'red'");
  gnuplot_plot_coordinates(h1, xp[0].data(), yp[0].data(), xp[0].size(), "inputs0");

  // Draw the green points.
  gnuplot_setstyle(h1, "points");
  gnuplot_append_style(h1, " pointtype 4 pointsize 2 linewidth 2 linecolor 'green'");
  gnuplot_plot_coordinates(h1, xp[1].data(), yp[1].data(), xp[1].size(), "inputs1");

  Vector xrange = line.xrange();
  Vector yrange = line.yrange();

  std::ostringstream xrange_ss;
  xrange_ss << "set xrange [" << xrange(0) << ':' << xrange(1) << ']';
  gnuplot_cmd(h1, xrange_ss.str().c_str());
  std::ostringstream yrange_ss;
  yrange_ss << "set yrange [" << yrange(0) << ':' << yrange(1) << ']';
  gnuplot_cmd(h1, yrange_ss.str().c_str());

  bool success = line.draw(h1, xrange(0), yrange(0), xrange(1), yrange(1));

  // Close plot window.
  if (success)
    std::cin.get();
  gnuplot_close(h1);
}

int main()
{
  std::random_device rd;
  RandomGenerator seed_generator(rd());
  std::uniform_int_distribution<int> seed_dist(0, 10000);
  int seed = seed_dist(seed_generator);
  std::cout << "seed = " << seed << std::endl;
  RandomGenerator generator(seed);

  Line line(generator);
  std::cout << "line = " << line << std::endl;

  std::vector<Vector> points = generate_points(generator, line, 5000);
  std::vector<Color> targets = generate_colors(generator, line, points);

  plot(line, points, targets);
}
