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

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Quadratic Bezier", 1200, 900);

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

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Relationship between P₀, P₁, θ, w, v and s", {},
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
    draw::TextStyle<> label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle solid_line_style{.line_color = color::black, .line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
    draw::ArcStyle arc_style{.line_color = color::blue, .line_width = 1.0};

//    for (int j = 0;; ++j)
    {
      // P₀, the point at t=0.
      Point P0(1.72, -6.74);
      auto plot_P0 = plot.create_point(second_layer, P0, point_style);
      auto P0_label = plot.create_text(second_layer, P0, "P₀", label_style);

      // P₁, the point at t=1.
      Point P1(2.28, -4.76);
      auto plot_P1 = plot.create_point(second_layer, P1, point_style);
      auto P1_label = plot.create_text(second_layer, P1, "P₁", label_style({.position = draw::centered_right_of}));

      // Draw a cirle around the midpoint of P₀P₁ with radius |P₀P₁|/2.
      Vector P0P1(P0, P1);
      Point P0P1_circle_center = P0 + 0.5 * P0P1;
      auto plot_P0P1_circle_center = plot.create_point(second_layer, P0P1_circle_center, point_style);
      double P0P1_circle_radius = 0.5 * P0P1.length();
      auto plot_P0P1_circle = plot.create_circle(second_layer, P0P1_circle_center, P0P1_circle_radius,
          line_style({.line_color = color::gray}));

      double phi = 2.52; // 0.01 * j;
      Point point_on_circle(P0P1_circle_center.x() + P0P1_circle_radius * cos(phi), P0P1_circle_center.y() + P0P1_circle_radius * sin(phi));
      auto plot_point_on_circle = plot.create_point(second_layer, point_on_circle, point_style);
      Vector P1R_P1(point_on_circle, P1);

      double sw = P1R_P1.length();
      double sign = (P1R_P1.dot(P0P1.rotate_90_degrees()) > 0.0) ? -1.0 : 1.0;
      double s_squared_times_one_minus_2v = sign * Vector(point_on_circle, P0).length();

      double s = 2.0;
      double w = sw / s;
      double v = 0.5 * (1 - s_squared_times_one_minus_2v / (s * s));

      // Determine the symmetry line direction.
      Direction const perpendicular_to_symmetry_line_dir = P1R_P1.direction();
      Direction const symmetry_line_dir = perpendicular_to_symmetry_line_dir.normal();
      double theta = perpendicular_to_symmetry_line_dir.as_angle();

      // Define the matrix M.
      double m00 = sw * std::cos(theta) + 2.0 * s * s * v * std::sin(theta);
      double m01 = -s * s * std::sin(theta);
      double m10 = sw * std::sin(theta) - 2.0 * s * s * v * std::cos(theta);
      double m11 = s * s * std::cos(theta);

      auto xt = [=](double t){ return P0.x() + t * (m00 + m01 * t); };
      auto yt = [=](double t){ return P0.y() + t * (m10 + m11 * t); };

      // Let t run from v-4 to v+4, and then plot P₀ + t M [1 t].
      std::vector<Point> curve_points;
      for (double t = v - 4.0; t <= v + 4.0; t += 0.01)
      {
        double x = xt(t);
        double y = yt(t);
        curve_points.emplace_back(x, y);
      }
      auto curve = plot.create_curve(second_layer, std::move(curve_points), curve_line_style);

      // V, the parabola vertex point resides at t=v.
      auto V = plot.create_point(second_layer, {xt(v), yt(v)}, point_style);
      label_style.position = draw::centered_left_of;
      auto V_label = plot.create_text(second_layer, V, "V", label_style({.position = draw::centered_below}));

      // Draw a vertical line from V up to y=-6.
      Point V6(V.x(), -6.0);
      LinePiece foo(V, V6);
      auto plot_foo = plot.create_line(second_layer, foo, line_style);

      // Point on symmetry line at distance 1 from V.
      Point V1 = V + symmetry_line_dir;
      auto plot_V1 = plot.create_point(second_layer, V1, point_style);
//      auto V1_label = plot.create_text(second_layer, V1, "V1", label_style);
      // Draw an arrow from V to V1.
      auto plot_V_V1 = plot.create_connector(second_layer, V, V1, solid_line_style({.line_color = color::purple}));

      // Draw the unit vectors X and Y.
      auto plot_X = plot.create_connector(second_layer, V, V + perpendicular_to_symmetry_line_dir, solid_line_style);
      auto X_label = plot.create_text(second_layer, V + 0.5 * perpendicular_to_symmetry_line_dir, "X", label_style);
      auto plot_Y = plot.create_connector(second_layer, V, V + symmetry_line_dir, solid_line_style);
      auto Y_label = plot.create_text(second_layer, V + 0.5 * symmetry_line_dir, "Y", label_style);

      // Go from V1 a distance w left and right.
      Point V1L = V1 + w * perpendicular_to_symmetry_line_dir;
      Point V1R = V1 - w * perpendicular_to_symmetry_line_dir;
      auto plot_V1L = plot.create_point(second_layer, V1L, point_style({.color_index = color_index2}));
      auto plot_V1R = plot.create_point(second_layer, V1R, point_style({.color_index = color_index2}));

      auto plot_P1R = plot.create_point(second_layer, P1 - sw * perpendicular_to_symmetry_line_dir, point_style({.color_index = color_index2}));
//      auto P1R_label = plot.create_text(second_layer, plot_P1R, "P1R", label_style({.position = draw::centered_below}));
      auto plot_P0_P1R = plot.create_connector(second_layer, P0, plot_P1R, line_style);
      auto P0_P1R_label = plot.create_text(second_layer, P0 + 0.5 * Vector(P0, plot_P1R), "s²(1-2v)",
          label_style({.font_size = 12, .offset = 5}));

      // Draw a line through P₀ and P₁.
      auto line_through_P0_and_P1 = plot.create_line(second_layer, P0, P1, solid_line_style);

      // Draw a line through P₁ perpendicular to the symmetry line of the parabola.
      auto line_through_p1_perpendicular_to_symmetry_line =
        plot.create_line(second_layer, P1, perpendicular_to_symmetry_line_dir, line_style);

      // Draw the symmetry line of the parabola, through V.
      auto symmetry_line_of_parabola = plot.create_line(second_layer, V, symmetry_line_dir, line_style);

      // Draw a line between V1 and V1L.
      auto line_at_one_from_V = plot.create_connector(second_layer,
          V1, V1L, Connector::open_arrow, Connector::open_arrow, line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));
      auto w_label = plot.create_text(second_layer, V1 + 0.5 * w * perpendicular_to_symmetry_line_dir,
          "w", label_style({.position = draw::centered_above, .font_size = 14, .offset = 5}));
      // Draw a line between P₁ and P1R.
      auto line_P1_P1R = plot.create_connector(second_layer,
          P1, plot_P1R, Connector::open_arrow, Connector::open_arrow, line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));
      auto sw_label2 = plot.create_text(second_layer, plot_P1R + 0.5 * sw * perpendicular_to_symmetry_line_dir,
          "sw", label_style({.position = draw::centered_below, .font_size = 12, .offset = 5}));

      // Draw an arc between P1_P1R and P1_P0.
      plot::Arc alpha_arc(P1, Direction{P1, plot_P1R}, Direction{P1, P0}, 0.1);
      plot.add_arc(second_layer, alpha_arc, arc_style({.line_color = color::blue}));
      auto alpha_label = plot.create_text(second_layer, P1 + 0.13 * alpha_arc.bisector_direction(),
          "α", label_style({.position = draw::centered, .font_size = 14}));

      plot::Arc theta_arc(V, Direction::up, symmetry_line_dir, 0.2);
      plot.add_arc(second_layer, theta_arc, arc_style);
      auto theta_label = plot.create_text(second_layer, V + 0.24 * theta_arc.bisector_direction(),
          "θ", label_style({.position = draw::centered, .font_size = 14}));

      window.set_send_expose_events(true);

      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      std::cout << "phi = " << phi << std::endl;
      std::cin.get();

      window.set_send_expose_events(false);
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
