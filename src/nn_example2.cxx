#include "DebugLayer.h"
#include <array>

int main()
{
  // We have two layers with each 2 inputs and 2 outputs (and a batch size of 1).
  std::array<Layer<2, 2, 1>, 2> nn;

  // Initialization.
  nn[0].set_weights({{
      { 0.15, 0.20, 0.35 },
      { 0.25, 0.30, 0.35 }
    }});
  nn[1].set_weights({{
      { 0.40, 0.45, 0.60 },
      { 0.50, 0.55, 0.60 }
    }});

  // Input to first layer.
  Matrix X{{2, 1}};
  X(0, 0) = 0.05;
  X(1, 0) = 0.10;

  // Target output.
  Matrix T{{2, 1}};
  T(0, 0) = 0.01;
  T(1, 0) = 0.99;

  for (int step = 0; step < 10000; ++step)
  {
    // Forward pass.
    nn[0].forward(X);
    nn[1].forward(nn[0].outputs());

    // Calculate loss (one for each batch sample).
    Vector L = mse(nn[1].outputs(), T);
    std::cout << "output = " << nn[1].outputs() << "; L = " << L << std::endl;
    //std::cout << "nn[0] = " << nn[0] << "; nn[1] = " << nn[1] << std::endl;

    auto xi = derivative_mse(nn[1].outputs(), T);
    xi = nn[1].back_propagate(0.5, xi);
    xi = nn[0].back_propagate(0.5, xi);
  }
}
