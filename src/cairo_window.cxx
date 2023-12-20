#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/draw/Rectangle.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Axes.h"
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

    // Open window, handle event loop and block until the window is closed.
    EventLoop event_loop = window.run();

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::gray COMMA_DEBUG_ONLY("background_layer"));

    // Draw something on the background layer.
    auto red_square = std::make_unique<draw::Rectangle>(Rectangle{50, 50, 350, 250}, color::red);
    background_layer->draw(red_square);

    // Create another layer.
    auto second_layer = window.create_layer<Layer>(Rectangle{100, 100, 600, 400} COMMA_DEBUG_ONLY("second_layer"));

    // Draw a line.
    draw::Line blue_line({350, 250, 100, 100}, color::blue);
    second_layer->draw(&blue_line);

    // Create and draw a green square.
    auto green_square = std::make_unique<draw::Rectangle>(Rectangle{350, 250, 100, 100}, color::green);
    second_layer->draw(green_square);

    // Create and draw plot axes.
    draw::Axes axes({49.5, 9.5, 701, 541}, color::black);
    background_layer->draw(&axes);

    // Wait for a key-press.
    std::cin.get();

    // Wait for the event loop thread to be finished (the window closed).
    event_loop.set_cleanly_terminated();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
