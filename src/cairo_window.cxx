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
    plot.set_xrange({-3, 14});
    plot.set_yrange({-12, 5});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    int color_index = color_pool.get_and_use_color();
    int filled_shape = 1;
    draw::PointStyle point_style(color_index, filled_shape);
    draw::TextStyle<> point_label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};

#if 0
    std::vector<Point> curve_points;
    for (double x = -3.0; x < 12.0; x += 0.1)
      curve_points.emplace_back(x, y(x));
#endif

    std::vector<Point> points;
    int l = 0;
#if 0
    char const* label[] = { "Q0", "P0", "P2", "Q2" };
    for (double xi = -3.0; xi <= 12.0; xi += 5.0)
    {
      double yi = y(xi);
      points.emplace_back(xi, yi);
      plot.add_point(second_layer, xi, yi, point_style);
      point_label_style.position = (l < 2) ? draw::centered_left_of : draw::centered_right_of;
      plot.add_text(second_layer, xi, yi, label[l], point_label_style);
      ++l;
    }
#endif

    draw::LineStyle curve_line_style{.line_width = 1.0};
#if 0
    plot.add_curve(second_layer, curve_points, curve_line_style);
#endif
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
#if 0
    plot.add_line(second_layer, points[0], points[2], line_style);
    plot.add_line(second_layer, points[1], points[3], line_style);

    line_style.dashes = {5.0, 5.0};
    line_style.line_width = 2.0;
    line_style.line_color = color::green;
    plot.add_line(second_layer, y(points[1].x() + 0.01) - y(points[1].x() - 0.01), -0.02, points[1], line_style);
    line_style.line_color = color::teal;
    plot.add_line(second_layer, y(points[2].x() + 0.01) - y(points[2].x() - 0.01), -0.02, points[2], line_style);
#endif

    //      ⎡   w^(2/3)          0⎤
    //  M = ⎣2s w^(-1/3)  w^(-2/3)⎦
    //
    // Let t run from -1 to 2, and then plot X₀ + t M [1 t].
    std::vector<Point> curve_points2;
    double P1x = 6.0;
    double P1y = -2.0;
    double P2x = 2.0;
    double P2y = -4.0;
    double w = 1.5;
    double s = 1.5;

    std::array<double, 1> angles = { /*0.0,*/ M_PI / 6 /*, M_PI / 2, M_PI*/ };
    for (double theta : angles)
    {
      double ct = std::cos(theta);
      double st = std::sin(theta);
      double w13 = std::pow(w, 1.0/3.0);
      double w23 = w13 * w13;

      double a = ct * w23 - 2.0 * st * s / w13;
      double b = -st / w23;
      double c = st * w23 + 2.0 * ct * s / w13;
      double d = ct / w23;
      auto xt = [=](double t){ return P1x + t * (a + b * t); };
      auto yt = [=](double t){ return P1y + t * (c + d * t); };
      for (double t = -6.0; t <= 2.0; t += 0.01)
      {
        double x = xt(t);
        double y = yt(t);
        curve_points2.emplace_back(x, y);
      }
      plot.add_curve(second_layer, curve_points2, curve_line_style);
      curve_points2.clear();

      // Vertex point.
      double vt = -s * w13;
      double vx = xt(vt);
      double vy = yt(vt);
      points.emplace_back(vx, vy);
      plot.add_point(second_layer, vx, vy, point_style);
      point_label_style.position = draw::centered_left_of;
      plot.add_text(second_layer, vx, vy, "V", point_label_style);
    }

    // P₁.
    points.emplace_back(P1x, P1y);
    plot.add_point(second_layer, P1x, P1y, point_style);
    point_label_style.position = draw::centered_left_of;
    plot.add_text(second_layer, P1x, P1y, "P₁", point_label_style);
    // P₂.
    points.emplace_back(P2x, P2y);
    plot.add_point(second_layer, P2x, P2y, point_style);
    point_label_style.position = draw::centered_left_of;
    plot.add_text(second_layer, P2x, P2y, "P₂", point_label_style);

#if 0
    // Unrotated curve y(x).
    auto Xy = [=](double x){ return P1y + (x - P1x) / a * (c + d * ((x - P1x) / a)); };

    // Draw vertical lines at Vx +/- w.
    plot.add_line(second_layer, 1.0, 0.0, Point{vx, Xy(vx)}, line_style);
    plot.add_line(second_layer, 1.0, 0.0, Point{vx + w, Xy(vx + w)}, line_style);
    plot.add_line(second_layer, 1.0, 0.0, Point{vx - w, Xy(vx - w)}, line_style);
    plot.add_line(second_layer, Point{vx, Xy(vx + w)}, Point{vx + w, Xy(vx + w)}, line_style);
    point_label_style.position = draw::centered_above;
    plot.add_text(second_layer, vx + 0.5 * w, Xy(vx + w), "w", point_label_style);
#endif

    // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
    // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
    EventLoop event_loop = window.run();

//    std::this_thread::sleep_for(std::chrono::seconds(1));
//    text.reset();

    event_loop.set_cleanly_terminated();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
