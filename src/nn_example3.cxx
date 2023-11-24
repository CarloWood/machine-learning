#include "DebugLayer.h"
#include <array>

int main()
{
  std::cout.precision(9);

  // We have one layer, 2 inputs, 1 output and a batch size of 9.
  Layer<2, 1, 9> nn;

  // Initialization.
  nn.set_weights({{
      { 0.30520469, 0.567137778, 0.418522716 }
    }});

  // Input to first layer.
  Matrix X{{2, 9}};
  std::array<FloatType, 9> x1 = { -4.40000725, -1.58161616, 1.4800446, 2.68721461, 5.12050915, 8.49706268, 10.5576601, 12.4890213, 16.4113674 };
  for (int s = 0; s < 9; ++s)
  {
    X(0, s) = s;
    X(1, s) = x1[s];
  }

  // Target output.
  Matrix T{{1, 9}};
  std::array<FloatType, 9> t = { 0, 1, 1, 0, 0, 1, 0, 0, 1 };
  for (int s = 0; s < 9; ++s)
    T(0, s) = t[s];

  std::cout << "nn = " << nn << std::endl;
  std::cout << "X = " << X << std::endl;
  std::cout << "T = " << T << std::endl;

  for (int step = 0; step < 5; ++step)
  {
    // Forward pass.
    nn.forward(X);

    // Calculate loss (one for each batch sample).
    Vector L = mse(nn.outputs(), T);
    std::cout << nn << std::endl;
    std::cout << "output = " << nn.outputs() << std::endl;
    std::cout << "L = " << L << std::endl;

    auto xi = derivative_mse(nn.outputs(), T);
    xi = nn.back_propagate(0.5, xi);
  }
}
