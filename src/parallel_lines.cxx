#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Vector.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include <thread>
#include <iostream>
#include "debug.h"

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
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Parallel lines", 1200, 900);

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

    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Parallel lines", {},
        "x", {},
        "y", {});
    plot.set_xrange({-2, 3});
    plot.set_yrange({-1, 4});
//    plot.set_xrange({-1.75, -0.5});
//    plot.set_yrange({-0.25, 0.75});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    int color_index = color_pool.get_and_use_color();
    int filled_shape = 10;
    draw::PointStyle point_style(color_index, filled_shape);
    draw::TextStyle<> point_label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0};
    draw::LineStyle dashed_line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
    draw::LineStyle D_line_style{.line_color = color::red, .line_width = 2.0, .dashes = {5.0, 5.0}};
    draw::LineStyle E_line_style{.line_color = color::blue, .line_width = 1.0};

    {
      // P₀
      Point P0(2.0, 2.0);

      double theta = M_PI / 7.0;
      double delta_theta = theta; // / 3.0;

      Direction const D0(theta);
      Direction const D1(theta + delta_theta);
      Line L0(D0, P0);
      Line vertical_line_m8(Direction::up, {-8.0 / 5, 0.0});
      Point I = L0.intersection_with(vertical_line_m8);
      Line L1(D1, I);
      Line vertical_line_6(Direction::up, {6.0 / 5, 0.0});
      Point P1 = L1.intersection_with(vertical_line_6);

      // P₀
      auto plot_P0 = plot.create_point(second_layer, P0, point_style);
      auto P0_label = plot.create_text(second_layer, P0, "P₀", point_label_style({.position = draw::centered_below}));
      // P₁
      auto plot_P1 = plot.create_point(second_layer, P1, point_style);
      auto P1_label = plot.create_text(second_layer, P1, "P₁", point_label_style);
      // I
//      auto plot_I = plot.create_point(second_layer, I, point_style);
//      auto I_label = plot.create_text(second_layer, I, "I", point_label_style);

      // L₀
      auto plot_L0 = plot.create_line(second_layer, L0, line_style);
      // L₁
      auto plot_L1 = plot.create_line(second_layer, L1, line_style);

      // Draw D₀
      Point D0_end_point = I + D0;
      plot::Connector plot_D0(I, D0_end_point);
      plot.add_connector(second_layer, plot_D0, D_line_style);
      auto D0_label = plot.create_text(second_layer, I + 0.75 * D0, "D0", point_label_style({.position = draw::centered_below}));
      // Draw D₁
      Point D1_end_point = I + D1;
      plot::Connector plot_D1(I, D1_end_point);
      plot.add_connector(second_layer, plot_D1, D_line_style);
      auto D1_label = plot.create_text(second_layer, I + 0.5 * D1, "D1", point_label_style({.position = draw::centered_above}));

      // The vector from P₀ to P₁.
      Vector P0P1(P0, P1);
      double len_P0P1 = P0P1.length();
      // Calculate R: the point where P₁ would be when rotating P₀P₁ onto L₀.
      Point R = P0 + len_P0P1 * D0;
      // Plot R with label.
      auto plot_R = plot.create_point(second_layer, R, point_style);
      auto label_R = plot.create_text(second_layer, R, "R", point_label_style({.position = draw::centered_below}));
      // The vector from P₁ to R.
      Vector P1R(P1, R);
      auto plot_P1_R = plot.create_connector(second_layer, P1, R, Connector::filled_arrow, E_line_style, color::blue);
      double len_P1R = P1R.length();
      double ratio = 0.5 * len_P1R / len_P0P1;
      static std::array<double, 6> coeffs = { 1.0, -1.0 / 2, -1.0 / 8, -1.0 / 16, -5.0 / 128, -7.0 / 256 };
      // double alpha = len_P1R * std::sqrt(1.0 - ratio * ratio);
      double alpha = len_P1R * taylor(ratio * ratio, coeffs);

      // Draw a line inbetween D₀ and D₁.
      auto bisector = plot.create_line(second_layer, P0, ((P1 - P0) + (R - P0)).direction(), line_style({.line_color = color::turquoise}));

      // Draw 𝛦.
      auto plot_E = plot.create_connector(second_layer, D0_end_point, D1_end_point, Connector::filled_arrow, E_line_style, color::blue);
      auto label_E = plot.create_text(second_layer, I + 0.5 * (D0 + D1), "E",
          point_label_style({.position = draw::centered_below}));
      // Draw a line through D1_end_point, perpendicular to L₀.
      auto plot_Ep = plot.create_line(second_layer, D1_end_point, D0.normal(), dashed_line_style({.line_color = color::blue}));
      // Same, through P₁.
      auto plot_P1p = plot.create_line(second_layer, P1, D0.normal(), dashed_line_style({.line_color = color::blue}));

      // Draw a line between P₀ and P₁.
      auto plot_line_P0_P1 = plot.create_line(second_layer, P0, P1, line_style);

      // Difference between D0 and D1.
      Vector E(D1.x() - D0.x(), D1.y() - D0.y());
      // The length of that vector.
      double len_E = E.length();
      Dout(dc::notice, "Expected |E| = " << (2 * std::sin(0.5 * delta_theta)));
      Dout(dc::notice, "Calculated |E| = " << len_E);
      double ratio2 = 0.5 * len_E;
      // double epsilon = len_E * std::sqrt(1.0 - ratio2 * ratio2);
      double epsilon = len_E * taylor(ratio2 * ratio2, coeffs);

      Dout(dc::notice, "Expected epsilon = " << std::tan(delta_theta));
      Dout(dc::notice, "Calculated epsilon = " << epsilon);

      Dout(dc::notice, "Expected alpha = " << (P0P1.dot(D0.normal())));
      Dout(dc::notice, "Calculated alpha = " << alpha);
      Point calc_I = P1 - (alpha / epsilon) * D1;
      auto plot_calc_I = plot.create_point(second_layer, calc_I, point_style);
      Dout(dc::notice, "calc_I = " << calc_I);

      // Wait until the window closed.
      event_loop.join();
    }
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
