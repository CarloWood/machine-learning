#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Vector.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <limits>
#include "debug.h"

long double debug_length_of(long double x, long double y)
{
  return std::sqrt(x * x + y * y);
}

template<size_t Size>
double taylor(double x, std::array<double, Size> const& coeffs)
{
  static_assert(Size > 0, "coeffs must have at least one element.");
  auto it = coeffs.rbegin();
  double result = *it++;
  while (it != coeffs.rend())
    result = x * result + *it++;
  return result;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Intersecting lines", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    bool joined = false;

//    std::random_device rd;
//    std::mt19937 engine(rd());
    for (int seed = 476; !joined; ++seed)
//    int seed = 10248;
    {
      Dout(dc::notice, "seed = " << seed);
      std::mt19937 engine(seed);

      // Determine the scale.
      std::uniform_real_distribution<double> scale_exp_dist(-6, 6);
      double scale_exp = scale_exp_dist(engine);
      double scale = std::pow(10.0, scale_exp);
      scale = 4.0;

      Dout(dc::notice, "scale = " << std::setprecision(std::numeric_limits<double>::max_digits10) << scale);

      // Generate the coordinates of the intersection point.
      std::uniform_real_distribution<double> scale_dist(-scale, scale);
      double Ix = scale_dist(engine);
      double Iy = scale_dist(engine);

      // Generate the coordinates of P0.
      double P0x = scale_dist(engine);
      double P0y = scale_dist(engine);

      // Switch to long double for high precision calculations.
      long double Ix_hp = Ix;
      long double Iy_hp = Iy;
      long double P0x_hp = P0x;
      long double P0y_hp = P0y;

      Dout(dc::notice, "hp: I = (" << std::setprecision(std::numeric_limits<double>::max_digits10) << Ix_hp << ", " << Iy_hp << ")");
      Dout(dc::notice, "hp: P0 = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) << P0x_hp << ", " << P0y_hp << ")");

      // Calculate D0 with high precision.
      long double I_P0_x_hp = P0x_hp - Ix_hp;
      long double I_P0_y_hp = P0y_hp - Iy_hp;
      long double I_P0_x_squared_hp = I_P0_x_hp * I_P0_x_hp;
      long double I_P0_y_squared_hp = I_P0_y_hp * I_P0_y_hp;
      long double len_I_P0_hp = std::sqrt(I_P0_x_squared_hp + I_P0_y_squared_hp);
      long double D0x_hp = I_P0_x_hp / len_I_P0_hp;
      long double D0y_hp = I_P0_y_hp / len_I_P0_hp;
      std::bernoulli_distribution distribution_50(0.5);
      if (distribution_50(engine))
      {
        D0x_hp = -D0x_hp;
        D0y_hp = -D0y_hp;
      }

      Dout(dc::notice, "hp: D0 = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) <<
          D0x_hp << ", " << D0y_hp << ")");
      Dout(dc::notice, "hp: |D0| = " << std::setprecision(std::numeric_limits<long double>::max_digits10) <<
          debug_length_of(D0x_hp, D0y_hp));

      // Determine Q with high precision (the projection of P1 onto L0).
      long double Qx_hp;
      long double Qy_hp;
      {
        double signed_len_Q = scale_dist(engine);
        long double signed_len_Q_hp = signed_len_Q;
        Qx_hp = Ix_hp + D0x_hp * signed_len_Q_hp;
        Qy_hp = Iy_hp + D0y_hp * signed_len_Q_hp;
      }

      Dout(dc::notice, "hp: Q = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) << Qx_hp << ", " << Qy_hp << ")");

      // Determine alpha; alpha is always positive.
      std::uniform_real_distribution<double> alpha_exp_dist(-8, -2);
//      std::uniform_real_distribution<double> alpha_exp_dist(-1, 0);
      double alpha_exp = alpha_exp_dist(engine);
      double alpha = scale * std::pow(10.0, alpha_exp);
//      alpha *= 30;

      // Calculate P1.
      long double P1x_hp;
      long double P1y_hp;
      {
        // Determine N0 with high precision.
        long double N0x_hp = -D0y_hp;
        long double N0y_hp = D0x_hp;

        // Calculate P1.
        long double alpha_hp = alpha;
        Dout(dc::notice, "hp: alpha = " << std::setprecision(std::numeric_limits<long double>::max_digits10) << alpha_hp);
        P1x_hp = Qx_hp + alpha_hp * N0x_hp;
        P1y_hp = Qy_hp + alpha_hp * N0y_hp;
      }
      Dout(dc::notice, "hp: P1 = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) << P1x_hp << ", " << P1y_hp << ")");

      long double P0_P1_x_hp = P1x_hp - P0x_hp;
      long double P0_P1_y_hp = P1y_hp - P0y_hp;
      Dout(dc::notice, "hp: P0P1 = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) << P0_P1_x_hp << ", " << P0_P1_y_hp << ")");
      long double len_P0_P1_hp = std::sqrt(P0_P1_x_hp * P0_P1_x_hp + P0_P1_y_hp * P0_P1_y_hp);
      Dout(dc::notice, "hp: |P0P1| = " << std::setprecision(std::numeric_limits<long double>::max_digits10) << len_P0_P1_hp);

      // R is the point on L0 where P1 ends up when we rotate it over the shortest distance towards L0.
      double Rx_hp;
      double Ry_hp;
      {
        if (P0_P1_x_hp * D0x_hp + P0_P1_y_hp * D0y_hp < 0.0)
        {
          Rx_hp = P0x_hp - len_P0_P1_hp * D0x_hp;
          Ry_hp = P0y_hp - len_P0_P1_hp * D0y_hp;
        }
        else
        {
          Rx_hp = P0x_hp + len_P0_P1_hp * D0x_hp;
          Ry_hp = P0y_hp + len_P0_P1_hp * D0y_hp;
        }
      }
      Dout(dc::notice, "hp: R = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) << Rx_hp << ", " << Ry_hp << ")");

      // Calculate D1 with high precision.
      long double D1x_hp;
      long double D1y_hp;
      {
        long double I_P1_x_hp = P1x_hp - Ix_hp;
        long double I_P1_y_hp = P1y_hp - Iy_hp;
        long double I_P1_x_squared_hp = I_P1_x_hp * I_P1_x_hp;
        long double I_P1_y_squared_hp = I_P1_y_hp * I_P1_y_hp;
        long double len_I_P1_hp = std::sqrt(I_P1_x_squared_hp + I_P1_y_squared_hp);
        D1x_hp = I_P1_x_hp / len_I_P1_hp;
        D1y_hp = I_P1_y_hp / len_I_P1_hp;
      }
      if (distribution_50(engine))
      {
        D1x_hp = -D1x_hp;
        D1y_hp = -D1y_hp;
      }

      Dout(dc::notice, "hp: D1 = (" << std::setprecision(std::numeric_limits<long double>::max_digits10) << D1x_hp << ", " << D1y_hp << ")");
      Dout(dc::notice, "hp: |D1| = " << std::setprecision(std::numeric_limits<long double>::max_digits10) << debug_length_of(D1x_hp, D1y_hp));

      if (D0x_hp * D1x_hp + D0y_hp * D1y_hp > 0)
      {
        Dout(dc::notice, "hp: |E| = " << std::setprecision(std::numeric_limits<long double>::max_digits10) <<
            debug_length_of(D1x_hp - D0x_hp, D1y_hp - D0y_hp));
      }
      else
      {
        Dout(dc::notice, "hp: |E| = " << std::setprecision(std::numeric_limits<long double>::max_digits10) <<
            debug_length_of(D1x_hp + D0x_hp, D1y_hp + D0y_hp));
      }

      // Define all input variables as double.
      double P1x = P1x_hp;
      double P1y = P1y_hp;
      double D0x = D0x_hp;
      double D0y = D0y_hp;
      double D1x = D1x_hp;
      double D1y = D1y_hp;

      // Create plot object.
      plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
          "Intersecting lines", {},
          "x", {},
          "y", {});
      // Set ranges.
      double D0x_ep = Ix + D0x;
      double D0y_ep = Iy + D0y;
      // Find min and max values.
      double min_x = std::min(std::min(std::min(std::min(std::min(P0x, P1x), D0x_ep), Ix), (double)Rx_hp), (double)Qx_hp);
      double max_x = std::max(std::max(std::max(std::max(std::max(P0x, P1x), D0x_ep), Ix), (double)Rx_hp), (double)Qx_hp);
      double min_y = std::min(std::min(std::min(std::min(std::min(P0y, P1y), D0y_ep), Iy), (double)Ry_hp), (double)Qy_hp);
      double max_y = std::max(std::max(std::max(std::max(std::max(P0y, P1y), D0y_ep), Iy), (double)Ry_hp), (double)Qy_hp);
      double half_range = std::max(max_x - min_x, max_y - min_y) * 0.55;
      double center_x = 0.5 * (min_x + max_x);
      double center_y = 0.5 * (min_y + max_y);
//      double center_x = Ix_hp;
//      double center_y = Iy_hp;
//      half_range = 0.1;
      Debug(dc::cairowindow.on());
      plot.set_xrange({center_x - half_range, center_x + half_range});
      plot.set_yrange({center_y - half_range, center_y + half_range});
      Debug(dc::cairowindow.off());
      plot.add_to(background_layer, true);

      // Define styles.
      utils::ColorPool<32> color_pool;
      int color_index = color_pool.get_and_use_color();
      int filled_shape = 10;
      draw::PointStyle point_style(color_index, filled_shape);
      draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0};
      draw::TextStyle<> point_label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
      draw::LineStyle dashed_line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
      draw::LineStyle D_line_style{.line_color = color::blue, .line_width = 2.0, .dashes = {5.0, 5.0}};
      draw::LineStyle E_line_style{.line_color = color::red, .line_width = 1.0};

      // Plot P0, P1, I and calc_I.
      Point P0{P0x, P0y};
      Point P1{P1x, P1y};
      Point I{Ix, Iy};
      Direction D0(Point{D0x, D0y});
      Direction D1(Point{D1x, D1y});
      Direction N1 = D1.normal();

//      if (D0.dot(N1) >= 0.01)
//        continue;
#if 1
      // Calculate I the standard way.

      // Take dot product of D0 with N1:
      double D0_dot_N1 = D0.dot(N1);
      Dout(dc::notice, "D0·N1 = " << D0_dot_N1);

      // Take dot product of P1P0 with N1:
      double P1P0_dot_N1 = (P1 - P0).dot(N1);
      Dout(dc::notice, "P1P0·N1 = " << P1P0_dot_N1);

      // Calculate lambda.
      double lambda = P1P0_dot_N1 / D0_dot_N1;
      Dout(dc::notice, "lambda = " << lambda);

      // Calculate intersection point.
      Point calc_I(P0x - lambda * D0x, P0y - lambda * D0y);
#endif

//      auto plot_circle = plot.create_circle(second_layer, P0, len_P0P1, line_style);

      Dout(dc::notice, "Calculated I = (" << std::setprecision(std::numeric_limits<double>::max_digits10) << calc_I << ")");

      // Create two lines.
      Line L0(D0, P0);
      Line L1(D1, P1);
      Point real_I = L0.intersection_with(L1);

      Dout(dc::notice, "Relative error in I = (" <<
          (100.0 * std::abs(calc_I.x() - Ix) / std::abs(Ix)) << "%, " <<
          (100.0 * std::abs(calc_I.y() - Iy) / std::abs(Iy)) << "%)");

      double scaled_abs_error_x = std::abs(calc_I.x() - Ix) / scale;
      double scaled_abs_error_y = std::abs(calc_I.y() - Iy) / scale;
      Dout(dc::notice, "Scaled absolute error in I = (" <<
          (100.0 * scaled_abs_error_x) << "%, " <<
          (100.0 * scaled_abs_error_y) << "%)");

      double real_scaled_abs_error_x = std::abs(real_I.x() - Ix) / scale;
      double real_scaled_abs_error_y = std::abs(real_I.y() - Iy) / scale;
      Dout(dc::notice, "Scaled absolute error in real_I = (" <<
          (100.0 * real_scaled_abs_error_x) << "%, " <<
          (100.0 * real_scaled_abs_error_y) << "%)");

      // P₀
      auto plot_P0 = plot.create_point(second_layer, P0, point_style);
      auto P0_label = plot.create_text(second_layer, P0, "P₀", point_label_style({.position = draw::centered_right_of}));
      // P₁
      auto plot_P1 = plot.create_point(second_layer, P1, point_style);
      auto P1_label = plot.create_text(second_layer, P1, "P₁", point_label_style);
      // I
      auto plot_I = plot.create_point(second_layer, I, point_style);
      auto I_label = plot.create_text(second_layer, I, "I", point_label_style);
      // calc_I
      auto plot_calc_I = plot.create_point(second_layer, calc_I, point_style({.color_index = 2}));
      auto calc_I_label = plot.create_text(second_layer, calc_I, "I", point_label_style({.position = draw::centered_right_of}));

      // Draw lines
      auto plot_line0 = plot.create_line(second_layer, P0, D0, line_style);
      auto plot_line1 = plot.create_line(second_layer, P1, D1, line_style);

      // Draw D₀
      Point D0_end_point = I + D0;
      plot::Connector plot_D0(I, D0_end_point);
      plot.add_connector(second_layer, plot_D0, D_line_style);
      auto D0_label = plot.create_text(second_layer, I + 0.5 * D0, "D0", point_label_style({.position = draw::centered_below}));
      // Draw D₁
      Point D1_end_point = I + D1;
      plot::Connector plot_D1(I, D1_end_point);
      plot.add_connector(second_layer, plot_D1, D_line_style);
      auto D1_label = plot.create_text(second_layer, I + 0.5 * D1, "D1", point_label_style({.position = draw::centered_above}));

      // Draw line through P0, parallel to N1.
      Line QP0(N1, P0);
      Point Q = QP0.intersection_with(plot_line1);
      auto plot_Q = plot.create_point(second_layer, Q, point_style({.color_index = 3}));
      auto plot_QP0 = plot.create_connector(second_layer, Q, P0, dashed_line_style);

      LinePiece ab(Q, P0);
      Dout(dc::notice, "|QP0| = " << std::setprecision(std::numeric_limits<double>::max_digits10) << ab.length());

      if (std::max(scaled_abs_error_x, scaled_abs_error_y) > 1e-7
          || std::max(real_scaled_abs_error_x, real_scaled_abs_error_y) > 1e-7)
      {
        // Wait until the window closed.
//        event_loop.join();
//        joined = true;
        std::cin.get();
        //break;
      }
    }
//    if (!joined)
//    {
//      event_loop.join();
//      joined = true;
//    }
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }
}
