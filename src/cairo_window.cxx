#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/draw/Rectangle.h"
#include "cairowindow/draw/Line.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
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
    Window window("My window", 800, 600);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

#if 0
    // Draw something on the background layer.
    auto red_square = std::make_unique<draw::Rectangle>(Rectangle{50, 50, 350, 250}, color::red);
    background_layer->draw(red_square);
#endif

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

#if 0
    // Draw a line.
    draw::Line blue_line({350, 250, 100, 100}, color::blue);
    second_layer->draw(&blue_line);

    // Create and draw a green square.
    auto green_square = std::make_unique<draw::Rectangle>(Rectangle{350, 250, 100, 100}, color::green);
    second_layer->draw(green_square);
#endif

    // Create and draw plot area.
//    draw::PlotArea plot_area({49.5, 9.5, 701, 541}, {.axes_line_width = 2.0});
//    plot_area.set_range(plot::x_axis, 0, 10);
//    plot_area.set_range(plot::y_axis, -10, 20);
//    background_layer->draw(&plot_area);

    // Add some text.
//    auto text = std::make_unique<draw::Text>("Hello world", draw::centered_below, 350, 350, draw::TextStyle{.font_size = 24.0});
//    second_layer->draw(text);

    plot::Plot plot(window.geometry(), {.grid = {.color = color::orange}},
        "This is the plot title", {},
        "x", {},
        "y", {});
    plot.set_xrange({0, 10});
    plot.set_yrange({-2.5, 2.5});
    plot.add_to(second_layer);

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
