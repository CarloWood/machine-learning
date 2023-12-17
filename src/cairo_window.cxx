#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
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

    // Open window, handle event loop and block until the window is closed.
    EventLoop event_loop = window.run();

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::gray);

    // Create another layer.
    auto second_layer = window.create_layer<Layer>(Rectangle{100, 100, 600, 400});

    // It doesn't really matter if you destroy the draw::Rectangle afterwards, at the moment.
    {
      auto green_square = second_layer->create_layer_region<draw::Rectangle>(Rectangle{350, 250, 100, 100}, color::green);
      green_square->draw();
    }

    // Draw something on the background layer.
    auto red_square = background_layer->create_layer_region<draw::Rectangle>(Rectangle{50, 50, 350, 250}, color::red);
    red_square->draw();

    // Draw a line.
    auto blue_line = second_layer->create_layer_region<draw::Line>(Rectangle{350, 250, 100, 100}, color::blue);
    blue_line->draw();

    // Wait for the event loop thread to be finished (the window closed).
    event_loop.set_cleanly_terminated();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
