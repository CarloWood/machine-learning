#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Connector.h"
#include "cairowindow/Vector.h"
#include "cairowindow/draw/Shape.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "cairowindow/draw/ArrowHead.h"
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
        "123y456", {});
    plot.set_xrange({-2, 2});
    plot.set_yrange({-2, 2});
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

    {
      // P₀
      Point P0(-6.0 / 5, 0.0);

      double theta = M_PI / 8.0;
      double delta_theta = theta / 3.0;
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
      auto P0_label = plot.create_text(second_layer, P0, "P₀", point_label_style);
      // P₁
      auto plot_P1 = plot.create_point(second_layer, P1, point_style);
      auto P1_label = plot.create_text(second_layer, P1, "P₁", point_label_style);
      // I
      auto plot_I = plot.create_point(second_layer, I, point_style);
      auto I_label = plot.create_text(second_layer, I, "I", point_label_style);

      // L₀
      auto plot_L0 = plot.create_line(second_layer, L0, line_style);
      // L₁
      auto plot_L1 = plot.create_line(second_layer, L1, line_style);

      // Draw D₀
      plot::Connector plot_D0(P0, P0 + D0);
      plot.add_connector(second_layer, plot_D0, D_line_style);
      auto D0_label = plot.create_text(second_layer, P0 + 0.5 * D0, "D0", point_label_style({.position = draw::centered_below}));
      // Draw D₁
      plot::Connector plot_D1(P0, P0 + D1);
      plot.add_connector(second_layer, plot_D1, D_line_style);
      auto D1_label = plot.create_text(second_layer, P0 + 0.5 * D1, "D1", point_label_style({.position = draw::centered_above}));

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
