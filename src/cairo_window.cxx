#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
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
    auto red_square = std::make_unique<draw::Shape>(Rectangle{50, 50, 350, 250},
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
//    auto text = std::make_unique<draw::Text>("Hello world", draw::centered_below, 350, 350, draw::TextStyle{.font_size = 24.0});
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
    int filled_shape = 1;
    draw::PointStyle point_style(color_index, filled_shape);
    draw::TextStyle<> point_label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};

    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};

//    for (int j = 0; j < 100; ++j)
//    for (int i = 0; i < 150; ++i)
    {
      // The point at t=0.
      double P0x = 2.0;
      double P0y = -6.0;

      // Curve characteristics.
      double w = 0.8;                     // "width"
      double s = 0.25; //-1.5 + i * 0.03;                    // "shift" (of P0 along the curve; if s=0 then P0 corresponds to the vertex point).
      double theta = 7 * M_PI / 6;            // Counter clock-wise rotation of the parabola in radians.

      // Define the matrix M.
      double m00 = w * std::cos(theta) + 2.0 * s * std::sin(theta);
      double m01 = -std::sin(theta);
      double m10 = w * std::sin(theta) - 2.0 * s * std::cos(theta);
      double m11 = std::cos(theta);

      auto xt = [=](double t){ return P0x + t * (m00 + m01 * t); };
      auto yt = [=](double t){ return P0y + t * (m10 + m11 * t); };

      // The point at t=1.
      double P1x = xt(1.0);
      double P1y = yt(1.0);

      // The distance between P₀ and P₁ is determined by w and s (see README.bezier):
      double P0P1 = std::sqrt(w * w + (2 * s - 1) * (2 * s - 1));

#if 1
      // Let t run from s-4 to s+4, and then plot P₀ + t M [1 t].
      std::vector<Point> curve_points2;
      for (double t = s - 4.0; t <= s + 4.0; t += 0.01)
      {
        double x = xt(t);
        double y = yt(t);
        curve_points2.emplace_back(x, y);
      }
      plot.add_curve(second_layer, curve_points2, curve_line_style);
      curve_points2.clear();
#endif

      // P₁.
      plot.add_point(second_layer, P1x, P1y, point_style);
      point_label_style.position = draw::centered_left_of;
      plot.add_text(second_layer, P1x, P1y, "P₁", point_label_style);

      // Vertex point.
      double vt = s;
      double vx = xt(vt);
      double vy = yt(vt);
      plot.add_point(second_layer, vx, vy, point_style);
      point_label_style.position = draw::centered_left_of;
      plot.add_text(second_layer, vx, vy, "V", point_label_style);

      // P₀.
      plot.add_point(second_layer, P0x, P0y, point_style);
      point_label_style.position = draw::centered_left_of;
      plot.add_text(second_layer, P0x, P0y, "P₀", point_label_style);
      // Draw a horizontal line through P₁.
      plot.add_line(second_layer, Point{P0x, P1y}, Point{P1x, P1y}, line_style);
      // Draw a line through P₀ and P₁.
      plot.add_line(second_layer, Point{P0x, P0y}, Point{P1x, P1y}, line_style);
      // Draw a line through P₁ perpendicular to the symmetry line of the parabola.
      plot.add_line(second_layer, -std::sin(theta), std::cos(theta), Point{P1x, P1y}, line_style);
      // Draw the symmetry line of the parabola, through V.
      plot.add_line(second_layer, std::cos(theta), std::sin(theta), Point{vx, vy}, line_style);
      // Draw a line perpendicular to the symmetry line of the parabola, at a distance of 1 from V.
      plot.add_line(second_layer, -std::sin(theta), std::cos(theta), Point{vx - std::sin(theta), vy + std::cos(theta)}, line_style);
      plot.add_point(second_layer, xt(5 * s), yt(5 * s), point_style);
      // Draw a vertical line through P₀.
      plot.add_line(second_layer, Point{P0x, P0y}, Point{P0x, P1y}, line_style);
//      plot.add_line(second_layer, 1.0, 0.0, Point{P1x, P1y}, line_style);
      // Draw a cirle around P₀ with radius P0P1.
//      plot.add_circle(second_layer, Point{P0x, P0y}, P0P1, line_style);

      //std::this_thread::sleep_for(std::chrono::milliseconds(100));
      std::cin.get();

      plot.remove_points();
      plot.remove_texts();
      plot.remove_lines();
      plot.remove_circles();
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
