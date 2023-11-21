#include <DebugTensor.h>
#include <iostream>
#include <format>

// Compile as: clang++ -std=c++20 -g -I. nn_example.cxx

double const epsilon = 0.000001;

int main()
{
  Function const a(sigmoid);
  Function const sqr(square);

  Vector input({3});

  input(0) = 0.05;
  input(1) = 0.1;
  input(2) = 1;

  Matrix wbih({2, 3});

  wbih(0, 0) = 0.15; wbih(0, 1) = 0.2; wbih(0, 2) = 0.35;
  wbih(1, 0) = 0.25; wbih(1, 1) = 0.3; wbih(1, 2) = 0.35;

  auto h01 = a(contract(wbih, input, 1, 0));
  Vector h({3});

  h(0) = h01(0);
  h(1) = h01(1);
  h(2) = 1;

  std::cout << "h = " << h << std::endl;

  Matrix wbho({2, 3});

  wbho(0, 0) = 0.40; wbho(0, 1) = 0.45; wbho(0, 2) = 0.6;
  wbho(1, 0) = 0.5;  wbho(1, 1) = 0.55; wbho(1, 2) = 0.6;

  auto o_a = a(contract(wbho, h, 1, 0));

  std::cout << "o_a = " << o_a << std::endl;

  Vector t({2});

  t(0) = 0.01;
  t(1) = 0.99;

  auto sd_a = sqr(t - o_a);
  double loss_a = 0.5 * (sd_a(0) + sd_a(1));

  std::cout << "loss_a = " << loss_a << std::endl;

#if 0
  // Emperically determine ∂L/∂zᵢ (our zᵢ is h here).

  for (int i = 0; i < 2; ++i)
  {
    Vector h_eps = h;
    h_eps(i) += epsilon;
    auto o_eps = a(contract(wbho, h_eps, 1, 0));
    auto sd_eps = sqr(t - o_eps);
    double loss_eps = 0.5 * (sd_eps(0) + sd_eps(1));

    std::cout << " ∂L/∂z_" << i << "=" << std::format("{: >10.7f}", (loss_eps - loss_a) / epsilon) << '\n';
  }
#endif

#if 0
  // Emperically determine ∂L/∂qᵢⱼ (our qᵢⱼ is wbho here).

  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      Matrix wbho_eps = wbho;
      wbho_eps(i, j) += epsilon;
      auto o_eps = a(contract(wbho_eps, h, 1, 0));
      auto sd_eps = sqr(t - o_eps);
      double loss_eps = 0.5 * (sd_eps(0) + sd_eps(1));

      std::cout << " ∂L/∂q_" << i << "," << j << "=" << std::format("{: >10.7f}", (loss_eps - loss_a) / epsilon);
    }
    std::cout << '\n';
  }
#endif

#if 0
  // Compose ∇_Q L (see comment in SigmoidDenseLayer in neural_network.cxx).

  Matrix A({1,2});
  A(0,0) = 0.74136507;
  A(0,1) = -0.21707153;

  Matrix B({2,2});
  B(0,0) = 0.18681560;
  B(1,1) = 0.17551005;

  Tensor<3> C({2,2,3});
  C(0,0,0) = 0.593269992;
  C(0,0,1) = 0.596884378;
  C(0,0,2) = 1;
  C(1,1,0) = 0.593269992;
  C(1,1,1) = 0.596884378;
  C(1,1,2) = 1;

  auto R = contract(contract(A, B, 1, 0), C, 1, 0);

  std::cout << R << std::endl;
#endif

#if 1
  // Emperically determine ∂L/∂wᵢⱼ (our wᵢⱼ is wbih here).

  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      Matrix wbih_eps = wbih;
      wbih_eps(i, j) += epsilon;

      auto h01_eps = a(contract(wbih_eps, input, 1, 0));
      Vector h_eps({3});

      h_eps(0) = h01_eps(0);
      h_eps(1) = h01_eps(1);
      h_eps(2) = 1;

      auto o_eps = a(contract(wbho, h_eps, 1, 0));

      auto sd_eps = sqr(t - o_eps);
      double loss_eps = 0.5 * (sd_eps(0) + sd_eps(1));

      std::cout << " ∂L/∂w_" << i << "," << j << "=" << std::format("{: >10.7f}", (loss_eps - loss_a) / epsilon);
    }
    std::cout << '\n';
  }
#endif

#if 0
  // Compose ∇_W L (see comment in SigmoidDenseLayer in neural_network.cxx).

  Matrix A({1,2});
  A(0,0) = 0.03635031;
  A(0,1) = 0.04137032;

  Matrix B({2,2});
  B(0,0) = 0.241300708;
  B(1,1) = 0.240613417;

  Tensor<3> C({2,2,3});
  C(0,0,0) = 0.05;
  C(0,0,1) = 0.1;
  C(0,0,2) = 1;
  C(1,1,0) = 0.05;
  C(1,1,1) = 0.1;
  C(1,1,2) = 1;

  auto R = contract(contract(A, B, 1, 0), C, 1, 0);

  std::cout << R << std::endl;
#endif

}
