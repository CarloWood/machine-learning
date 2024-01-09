#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Vector.h"
#include "cairowindow/draw/Shape.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include <thread>
#include <iostream>
#include "debug.h"

double y(double x)
{
  return 2.5 - (x - 5) * (x - 5) / 5;
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
    Window window("My window", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

#if 0
    // Draw something on the background layer.
    auto red_square = std::make_shared<draw::Shape>(Rectangle{50, 50, 350, 250},
        draw::ShapeStyle{.line_color = color::red, .shape = draw::rectangle});
    background_layer->draw(red_square);
#endif

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

#if 0
    // Draw a line.
    draw::Line blue_line({350, 250, 100, 100}, draw::LineStyle{.line_color = color::blue, .line_width = 1.0});
    second_layer->draw(&blue_line);
#endif

    // Create and draw plot area.
//    draw::PlotArea plot_area({49.5, 9.5, 701, 541}, {.axes_line_width = 2.0});
//    plot_area.set_range(plot::x_axis, 0, 10);
//    plot_area.set_range(plot::y_axis, -10, 20);
//    background_layer->draw(&plot_area);

    // Add some text.
//    auto text = std::make_shared<draw::Text>("Hello world", draw::centered_below, 350, 350, draw::TextStyle{.font_size = 24.0});
//    second_layer->draw(text);

    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Constructing a Bezier curve", {},
        "x", {},
        "y", {});
    plot.set_xrange({0, 4});
    plot.set_yrange({-8, -4});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    int color_index = color_pool.get_and_use_color();
    int color_index2 = 3; //color_pool.get_and_use_color();
    Dout(dc::notice, "color_index2 = " << color_index2);
    int filled_shape = 1;
    draw::PointStyle point_style(color_index, filled_shape);
    draw::TextStyle<> point_label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle solid_line_style{.line_color = color::black, .line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};

//    for (int j = 0; j < 100; ++j)
    for (int i = 50; i < 150; ++i)
    {
      // P₀, the point at t=0.
      plot::Point P0(2.0, -6.0);

      // Curve characteristics.
      double w = 0.8;                     // "width"
      double s = /*0.25;*/ -1.5 + i * 0.03;                    // "shift" (of P0 along the curve; if s=0 then P0 corresponds to the vertex point).
      double theta = M_PI / 6;            // Counter clock-wise rotation of the parabola in radians.
      Direction const symmetry_line_dir(0.5 * M_PI + theta);
      Direction const perpendicular_to_symmetry_line_dir(theta);

      // Define the matrix M.
      double m00 = w * std::cos(theta) + 2.0 * s * std::sin(theta);
      double m01 = -std::sin(theta);
      double m10 = w * std::sin(theta) - 2.0 * s * std::cos(theta);
      double m11 = std::cos(theta);

      auto xt = [=](double t){ return P0.x() + t * (m00 + m01 * t); };
      auto yt = [=](double t){ return P0.y() + t * (m10 + m11 * t); };

      // Let t run from s-4 to s+4, and then plot P₀ + t M [1 t].
      std::vector<Point> curve_points;
      for (double t = s - 4.0; t <= s + 4.0; t += 0.01)
      {
        double x = xt(t);
        double y = yt(t);
        curve_points.emplace_back(x, y);
      }
      auto curve = plot.create_curve(second_layer, std::move(curve_points), curve_line_style);

      // P₁, the point at t=1.
      auto P1 = plot.create_point(second_layer, {xt(1.0), yt(1.0)}, point_style);
      point_label_style.position = draw::centered_left_of;
      auto P1_label = plot.create_text(second_layer, P1, "P₁", point_label_style({.position = draw::centered_right_of}));

      // V, the parabola vertex point resides at t=s.
      auto V = plot.create_point(second_layer, {xt(s), yt(s)}, point_style);
      point_label_style.position = draw::centered_left_of;
      auto V_label = plot.create_text(second_layer, V, "V", point_label_style({.position = draw::centered_below}));

      // Point on symmetry line at distance 1 from V.
      Point V1 = V + symmetry_line_dir;
      auto plot_V1 = plot.create_point(second_layer, V1, point_style);
      auto V1_label = plot.create_text(second_layer, V1, "V1", point_label_style);

      // Go from V1 a distance w left and right.
      Point V1L = V1 + w * perpendicular_to_symmetry_line_dir;
      Point V1R = V1 - w * perpendicular_to_symmetry_line_dir;
      auto plot_V1L = plot.create_point(second_layer, V1L, point_style({.color_index = color_index2}));
      auto plot_V1R = plot.create_point(second_layer, V1R, point_style({.color_index = color_index2}));

      auto plot_P1R = plot.create_point(second_layer, P1 - w * perpendicular_to_symmetry_line_dir, point_style({.color_index = color_index2}));
      auto P1R_label = plot.create_text(second_layer, plot_P1R, "P1R", point_label_style({.position = draw::centered_below}));
      auto plot_P0_P1R = plot.create_connector(second_layer, P0, plot_P1R, line_style);
      auto P0_P1R_label = plot.create_text(second_layer, P0 + 0.5 * Vector(P0, plot_P1R), "1-2s", point_label_style({.font_size = 12, .offset = 5}));

      // Helper point.
//      auto H = plot.create_point(second_layer, {P0.x(), P1.y()}, point_style);
//      auto H_label = plot.create_text(second_layer, H, "H", point_label_style);

      // P₀.
      plot.add_point(second_layer, P0, point_style);
      point_label_style.position = draw::centered_left_of;
      auto P0_label = plot.create_text(second_layer, P0, "P₀", point_label_style);

      // Draw a cirle around the midpoint of P₀P₁ with radius |P₀P₁|/2.
      Vector P0P1(P0, P1);
      auto circle = plot.create_circle(second_layer, P0 + 0.5 * P0P1, P0P1.length() / 2, line_style);

      // Draw a line through P₁ and H (horizontal because H has the same y coordinate as P₁).
//      auto horizontal_line_through_P1_and_H = plot.create_line(second_layer, P1, H, line_style);

      // Draw a line through P₀ and P₁.
      auto line_through_P0_and_P1 = plot.create_line(second_layer, P0, P1, solid_line_style);

      // Draw a line through P₁ perpendicular to the symmetry line of the parabola.
      auto line_through_p1_perpendicular_to_symmetry_line =
        plot.create_line(second_layer, P1, perpendicular_to_symmetry_line_dir, line_style);

      // Draw the symmetry line of the parabola, through V.
      auto symmetry_line_of_parabola =
        plot.create_line(second_layer, V, symmetry_line_dir, line_style);

      // Draw a line between V1 and V1L.
      auto line_at_one_from_V = plot.create_connector(second_layer,
          V1, V1L, Connector::open_arrow, Connector::open_arrow, line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));
      auto w_label = plot.create_text(second_layer, V1 + 0.5 * w * perpendicular_to_symmetry_line_dir,
          "w", point_label_style({.position = draw::centered_above}));
      // Draw a line between P₁ and P1R.
      auto line_P1_P1R = plot.create_connector(second_layer,
          P1, plot_P1R, Connector::open_arrow, Connector::open_arrow, line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));
      auto w_label2 = plot.create_text(second_layer, plot_P1R + 0.5 * w * perpendicular_to_symmetry_line_dir,
          "w", point_label_style({.position = draw::centered_above}));

      // Draw a line through P₀ and H (vertical because H has the same x coordinate as P₀).
//      auto vertical_line_through_P0 = plot.create_line(second_layer, P0, H, line_style);

      //std::this_thread::sleep_for(std::chrono::milliseconds(100));
      std::cin.get();
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
