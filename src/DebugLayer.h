#pragma once

#include "DebugTensor.h"

template<int n_inputs, int n_outputs, int batch_size>
class Layer
{
 public:
  Matrix inputs_{{n_inputs, batch_size}};
  Matrix outputs_{{n_outputs, batch_size}};

  Matrix weights_{{ n_outputs, n_inputs + 1 }};
  Function activation{sigmoid};

 public:
  void set_weights(std::array<std::array<FloatType, n_inputs + 1>, n_outputs> weights)
  {
    for (int i = 0; i < n_outputs; ++i)
      for (int j = 0; j < n_inputs + 1; ++j)
        weights_(i, j) = weights[i][j];
  }

  void forward(Matrix inputs)
  {
    Vector ones{{batch_size}, 1};
    inputs_ = inputs.append(0, ones);
    outputs_ = activation(contract(weights_, inputs_, 1, 0));
  }

  Matrix const& inputs() const { return inputs_; }
  Matrix const& outputs() const { return outputs_; }

  Matrix back_propagate(FloatType alpha, Matrix const& xi_l_plus_1)
  {
    Matrix delta{{n_outputs, batch_size}};
    Matrix xi_l{{n_outputs, batch_size}};
    for (int s = 0; s < batch_size; ++s)
    {
      // δₗᵢ = ξ₍ₗ₊₁₎ᵢ zᵢ(1-zᵢ)
      for (int i = 0; i < n_outputs; ++i)
        delta(i, s) = xi_l_plus_1(i, s) * outputs_(i, s) * (1.0 - outputs_(i, s));
      // ξₗᵢ = \sum_{k=0}^{Mₗ-1}(δₗₖ wₖᵢ)
      for (int i = 0; i < n_outputs; ++i)
      {
        FloatType sum = 0;
        for (int k = 0; k < n_outputs; ++k)
          sum += delta(k, s) * weights_(k, i);
        xi_l(i, s) = sum;
      }
    }
    for (int s = 0; s < batch_size; ++s)
    {
      // wᵢⱼ' = wᵢⱼ - α δₗᵢ xⱼ
      for (int i = 0; i < n_outputs; ++i)
        for (int j = 0; j < n_inputs + 1; ++j)
          weights_(i, j) -= (alpha / batch_size) * delta(i, s) * inputs_(j, s);
    }

    return xi_l;
  }

  void print_on(std::ostream& os) const
  {
    os << "{weights_:" << weights_ << '}';
  }

  friend std::ostream& operator<<(std::ostream& os, Layer const& layer)
  {
    layer.print_on(os);
    return os;
  }
};
