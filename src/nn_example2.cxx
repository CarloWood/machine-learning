#include "DebugTensor.h"

class Layer
{
 public:
  Vector inputs_{{2}};
  Vector outputs_{{2}};

  Matrix weights_{{ 2, 3 }};
  Function activation{sigmoid};

 public:
  void set_weights(std::array<std::array<double, 3>, 2> weights)        // A [2x3] matrix.
  {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 3; ++j)
        weights_(i, j) = weights[i][j];
  }

  void forward(Vector inputs)
  {
    inputs_ = inputs;
    outputs_ = activation(contract(weights_, inputs.append(1), 1, 0));
  }

  Vector const& outputs() const { return outputs_; }

  Vector back_propagate(double alpha, Vector const& xi_l_plus_1)
  {
    // δₗᵢ = ξ₍ₗ₊₁₎ᵢ zᵢ(1-zᵢ)
    Vector delta{{2}};
    for (int i = 0; i < 2; ++i)
      delta(i) = xi_l_plus_1(i) * outputs_(i) * (1.0 - outputs_(i));
    // ξₗᵢ = \sum_{k=0}^{Mₗ-1}(δₗₖ wₖᵢ)
    Vector xi_l{{2}};
    for (int i = 0; i < 2; ++i)
    {
      double sum = 0;
      for (int k = 0; k < 2; ++k)
        sum += delta(k) * weights_(k, i);
      xi_l(i) = sum;
    }
    // wₖᵢ' = wₖᵢ - α δₗᵢ xⱼ

    return xi_l;
  }

  void print_on(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, Layer const& layer)
  {
    layer.print_on(os);
    return os;
  }
};

void Layer::print_on(std::ostream& os) const
{
  os << "{weights_:" << weights_ << '}';
}

int main()
{
  // We have two layers.
  std::array<Layer, 2> nn;

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
  Vector X{{2}};
  X(0) = 0.05;
  X(1) = 0.10;

  // Target output.
  Vector T{{2}};
  T(0) = 0.01;
  T(1) = 0.99;

  for (int step = 0; step < 100; ++step)
  {
    // Forward pass.
    nn[0].forward(X);
    nn[1].forward(nn[0].outputs());

    // Calculate loss.
    double L = mse(nn[1].outputs(), T);
    std::cout << "output = " << nn[1].outputs() << "; L = " << L << std::endl;

    auto xi = derivative_mse(nn[1].outputs(), T);

    xi = nn[1].back_propagate(xi);
    xi = nn[0].back_propagate(xi);
  }
}
