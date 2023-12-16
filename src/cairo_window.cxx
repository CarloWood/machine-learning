#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
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

    // Create a new layer with a white background.
    auto background_layer = window.create_background_layer<Layer>(Color{1, 1, 1});

    // Create another layer.
    auto second_layer = window.create_layer<Layer>(Rectangle{100, 100, 600, 400});

    auto green_square = second_layer->create_layer_region<LayerRegion>();
    green_square->draw([](cairo_t* cr) -> Rectangle {
      cairo_set_source_rgb(cr, 0, 255, 0); // Green color for drawing.
      cairo_set_line_width(cr, 2);
      cairo_move_to(cr, 350, 250);
      cairo_line_to(cr, 350, 350);
      cairo_line_to(cr, 450, 350);
      cairo_line_to(cr, 450, 250);
      cairo_line_to(cr, 350, 250);
      cairo_stroke(cr);
      return {349, 249, 102, 102};
    });

    // Draw something on the background layer.
    auto red_square = background_layer->create_layer_region<LayerRegion>();
    red_square->draw([](cairo_t* cr) -> Rectangle {
      cairo_set_source_rgb(cr, 255, 0, 0); // Red color for drawing.
      cairo_set_line_width(cr, 2);
      cairo_move_to(cr, 50, 50);
      cairo_line_to(cr, 400, 50);
      cairo_line_to(cr, 400, 300);
      cairo_line_to(cr, 50, 300);
      cairo_line_to(cr, 50, 50);
      cairo_stroke(cr);
      return {49, 49, 352, 252};    // This area includes the line width.
    });

    // Wait for the event loop thread to be finished (the window closed).
    event_loop.set_cleanly_terminated();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
