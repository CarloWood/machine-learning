#include "gnuplot_i.h"
#include "DebugShow.h"
#include "dump_graph.h"
#include <fstream>
#include <random>
#include <format>

using namespace tensorflow;
using namespace tensorflow::ops;

constexpr int size_of_batch = 100;

float x1(float x0)
{
  return 2.5f * x0 - 4.3f;
}

bool AreTensorsEqual(tensorflow::Tensor const& tensor_a, tensorflow::Tensor const& tensor_b)
{
  // Check if shapes are the same.
  if (!tensor_a.shape().IsSameSize(tensor_b.shape()))
    return false;

  // Assuming tensors are of type float
  auto flat_a = tensor_a.flat<float>();
  auto flat_b = tensor_b.flat<float>();

  for (int i = 0; i < flat_a.size(); ++i)
    if (flat_a(i) != flat_b(i))
      return false;

  return true;
}

void AssertEqual(Scope const& scope, tensorflow::Tensor const& tensor_a, tensorflow::Tensor const& tensor_b)
{
  if (!AreTensorsEqual(tensor_a, tensor_b))
  {
    std::cerr << "Error: Tensors are not equal, aborting!";
    Abort(scope, Abort::Attrs().ErrorMsg("Error: Tensors are not equal, aborting!"));
  }
}

void draw_points(gnuplot_ctrl* h1, std::array<std::vector<double>, 2> const& xp, std::array<std::vector<double>, 2> const& yp)
{
  gnuplot_append_style(h1, " pointtype 7 pointsize 2 linewidth 2 linecolor 'red'");
  gnuplot_plot_coordinates(h1, xp[0].data(), yp[0].data(), xp[0].size(), "inputs0");
  gnuplot_setstyle(h1, "points");
  gnuplot_append_style(h1, " pointtype 4 pointsize 2 linewidth 2 linecolor 'green'");
  gnuplot_plot_coordinates(h1, xp[1].data(), yp[1].data(), xp[1].size(), "inputs1");
  gnuplot_cmd(h1, "unset key");
  gnuplot_cmd(h1, "set yrange [-15:20]");
  gnuplot_setstyle(h1, "points");
}

int main()
{
  std::cout.precision(9);

  std::random_device rd;
  std::default_random_engine seed_generator(rd());

  std::uniform_int_distribution<int> seed_distribution(0, 1000000000);
  int seed1 = seed_distribution(seed_generator);
  int seed2 = seed_distribution(seed_generator);

//  seed1 = 136637526; seed2 = 806889242;

  std::cout << "seed1 = " << seed1 << "; seed2 = " << seed2 << '\n';

  // GNU plot handle.
  gnuplot_ctrl* h1 = gnuplot_init();

  // Define the scope.
  Scope scope = Scope::NewRootScope().ExitOnError();
  DebugShow debug_show(scope);
  std::vector<Operation> initialization_operations;

  //         .---.
  // x₀ ---> |   | ---> [tanh activation] ---> { >0 --> label = 1.
  // x₁ ---> |   |                             { <0 --> label = 0.
  //         '---'
  //
  // Notation: row: [...], column: <...>, depth: {...}. All indices are in the order row, column, depth.
  // X is the input column-vector <x₀ x₁ 1>.
  // W is the weights/bias row-vector [w₀ w₁ b].
  //
  // Pre-activation output = WX = w₀x₀ + w₁x₁ + b.
  // Since the activation function tanh leaves the sign unchanged, we want the pre-activation output
  // to be zero when a point is exactly on the line, positive when it is above the line, and negative
  // when it is below the line. This means that w₁ must be positive.
  //
  // Using the above function x1(x0): x₁ = 2.5 x₀ - 4.3,
  // we have: w₀x₀ + w₁(2.5 x₀ - 4.3) + b = 0 -->
  // (w₀ + 2.5 w₁)x₀ - 4.3⋅w₁ + b = 0 for any x₀ -->
  // w₀ + 2.5⋅w₁ = 0 and -4.3⋅w₁ + b = 0 -->
  // w₀ = -2.5⋅w₁, and b = 4.3⋅w₁.
  // For example, if w₁=1 > 0, then w₀=-2.5 and b=4.3.
  //
  // This proves that a single perceptron with two inputs is able to learn
  // how to separate 2D points from laying above or below a line.

  // Create a neural network layer existing of a single perceptron with two inputs and one output.
  constexpr int input_units = 2;
  constexpr int output_units = 1;
  constexpr int UnknownRank = -1;       // Dynamic batch size.
  constexpr float alpha = 0.01;

  //---------------------------------------------------------------------------
  // Create graph.

  // Define the input placeholder with shape [2, UnknownRank], where 'UnknownRank' means the first dimension is dynamic.
  auto inputs = Placeholder(scope.WithOpName("Inputs"), DT_FLOAT, Placeholder::Shape({input_units, UnknownRank}));
  debug_show("inputs", inputs);

  // Define the labels placeholder with shape [UnknownRank], where 'UnknownRank' means the dimension is dynamic.
  auto labels = Placeholder(scope.WithOpName("Labels"), DT_FLOAT, Placeholder::Shape({UnknownRank}));
  debug_show("labels", labels);

  // Define a weights_bias matrix.
  auto weights_bias = Variable(scope.WithOpName("Weights"), {output_units, input_units + 1}, DT_FLOAT);

  // Define a persistent inputs_1 Variable.
  auto inputs_1 = Variable(scope.WithOpName("Inputs1"), {input_units + 1, size_of_batch}, DT_FLOAT);

  // Define a persistent variable for the batch_size, so we don't have to feed inputs
  // every epoch just to get the batch size from it.
  auto batch_size = Variable(scope.WithOpName("BatchSize"), {1}, DT_INT32);

  //---------------------
  // Initialization

  // Initialization of the weights and bias.
  initialization_operations.push_back(Assign(scope.WithOpName("assignW"), weights_bias,
        RandomNormal(scope.WithOpName("Rand"), {output_units, input_units + 1}, DT_FLOAT, RandomNormal::Attrs().Seed(seed1))).operation);
  debug_show("weights_bias", weights_bias);

  // Extract and save the batch size.
  auto batch_size_tensor = Slice(scope.WithOpName("BatchSizeTensor"), Shape(scope, inputs), {1}, {1});
  debug_show("batch_size_tensor", batch_size_tensor);

  // Initialize the variable batch_size with the actual batch size from inputs.
  Operation batch_size_assign_op = Assign(scope, batch_size, batch_size_tensor).operation;
  initialization_operations.push_back(batch_size_assign_op);
  debug_show(scope.WithControlDependencies({batch_size_assign_op}), "batch_size", batch_size);

  // Append a 1 to all inputs.
  auto inputs_1_tensor = Concat(scope.WithOpName("Inputs1Tensor"), {Input(inputs),
      ExpandDims(scope, Fill(scope, batch_size_tensor, 1.0f), 0)}, 0);
  debug_show("inputs_1_tensor", inputs_1_tensor);

  // Initialization of inputs_1 with the value corresponding to inputs.
  initialization_operations.push_back(Assign(scope, inputs_1, inputs_1_tensor).operation);
  debug_show("inputs_1", inputs_1);

  // Divide alpha by batch_size for use during backpropagation.
  auto learning_rate = Div(scope.WithOpName("Rate"), alpha, Cast(scope, batch_size, DT_FLOAT));
  debug_show("learning_rate", learning_rate);

  //---------------------
  // Forward propagation.

  // Because w₀₁ must be positive in order to uniquely define what is "above" and "below" the line,
  // change the sign of W if w₀₁ is negative.
  // First get W_01; Reshape returns a rank-2 1x1 matrix, reshape it to a rank-1 vector (can't reshape to a scalar directly).
  auto W_01 = Reshape(scope, Slice(scope, weights_bias, {0, 1}, {1, 1}), {1});
  // Create a rank-1 vector (of size 1) out of that, that contains true if W_01 is less than 0.
  auto W_01_is_negative_rank1 = Less(scope, W_01, {0.0f});
  // Now reduce this to a scalar, because Switch requires a scalar boolean as pred.
  auto W_01_is_negative = ReduceAny(scope, W_01_is_negative_rank1, 0);
  // Depending on the sign continue with either switch_op.output_false or switch_op.output_true.
  auto switch_op = Switch(scope, weights_bias, W_01_is_negative);
  // If switch_op.output_true was set, then W_01 is negative: negate all of weights_bias.
  Output negated_weights_bias = Negate(scope, switch_op.output_true);
  // Move the result into abs_weights_bias.
  auto abs_weights_bias = Merge(scope, { switch_op.output_false, negated_weights_bias }).output;
  // Update weights_bias with the possibly negated result.
  auto update_weights_bias_op2 = Assign(scope, weights_bias, abs_weights_bias, Assign::Attrs().UseLocking(true)).operation;

  // Perform `weights_bias * inputs_1` and apply activation function.
  auto z = MatMul(scope.WithControlDependencies({update_weights_bias_op2}), weights_bias, inputs_1);

  debug_show("z", z);
  auto outputs = Sigmoid(scope, z);
  debug_show("outputs", outputs);

  // Calculate the difference between the outputs and the targets.
  auto residual = Subtract(scope.WithOpName("Residual"), outputs, labels);     // Shape: [output_units x batch_size]

#if 0
  // Define a loss function, lets use Mean Squared Error (MSE) (totally random).
  // The loss function is not really used, except for printing.
  auto loss = ReduceMean(scope.WithOpName("Loss"), Square(scope, residual), {0});
#else
  // Use binary cross-entropy (BCE) loss function (using t = targets = labels):
  // loss = -(t log(outputs) + (1 - t) log(1 - outputs))
//  auto loss = ReduceMean(scope.WithOpName("Loss"), Negate(scope, Add(scope, Xlogy(scope, labels, outputs), Xlogy(scope, Sub(scope, 1.f, labels), Sub(scope, 1.f, outputs)))), {0});

  // Combine the binary cross-entropy loss (BCE) function with the sigmoid activation function (go straight from z to loss)
  // for greater numerical stability (using that outputs = sigmoid(z).
  //
  // loss = -t log(sigmoid(z)) - (1 - t) log(1 - sigmoid(z)) =
  //      = -t log(1 / (1 + exp(-z))) - (1 - t) log(1 - 1 / (1 + exp(-z))) =
  //      = t log(1 + exp(-z)) + (1 - t) log(1 + exp(z)) =
  //      = max(z, 0) - t * z + log(1 + exp(-abs(z)))
  auto loss = Add(scope,
      Sub(scope, Maximum(scope, z, 0.f),
                 Multiply(scope, labels, z)),
      Log(scope, Add(scope, 1.f, Exp(scope, Negate(scope, Abs(scope, z))))));
#endif
  debug_show("loss", loss);

  //---------------------
  // Backward propagation.

  // Calculate the initial xi value. This is ∂L/∂zᵢ, where L is the loss and Z the output (weights_bias * inputs_1).
  // So really this is 2/M (output - labels), but since M = output_units = 1 that is just 2 * residual.
  //
  // As described in README.back_propagation:
  //
  //    ξ₍ₗ₊₁₎ᵢ = 2/Mₗ (zᵢ - tᵢ)
  //
  // where i runs over the outputs (zᵢ) of the last layer (0 <= i < output_units (Mₗ)).
  //
  // Adding an extra dimension for the batch size, and leaving away the l index
  // that stands for layer, we have:
  //
  //    ξᵢₛ = 2 residualᵢₛ
  //
  // where s runs over all the samples in the batch.

  // Derivative of the loss function.
#if 0
  auto xi_1 = Multiply(scope.WithOpName("Xi_1"), 2.f, residual);   // Shape: [output_units x batch_size]
  debug_show("xi_1", xi_1);

  // Calculate delta as xi_1 times the derivative of the activation function.
  // For the sigmoid function we have: sigmoid' = outputs * (1 - outputs).
  //
  // As described in README.back_propagation:
  //
  //    δₗᵢ = ξ₍ₗ₊₁₎ᵢ zᵢ(1-zᵢ)
  //
  // where i runs over the outputs (zᵢ) of the current layer (0 <= i < output_units).
  //
  // Again adding an extra dimension for the batch size and leaving l away, we have:
  //
  //    δᵢₛ = ξᵢₛ zᵢₛ(1-zᵢₛ)
  //
  // where s runs over all the samples in the batch.
  //
  auto derivative_activation = Multiply(scope, outputs, Sub(scope, 1.f, outputs));
  auto delta = Multiply(scope.WithOpName("Delta"), xi_1, derivative_activation);      // Shape: [output_units x batch_size]
#else
  //auto xi_1 = Negate(scope.WithOpName("Xi_1"), Sub(scope, Xdivy(scope, labels, outputs), Xdivy(scope, Sub(scope, 1.f, labels), Sub(scope, 1.f, outputs))));
  // Calculate delta as Derivative of the combined activation + loss function.
  auto delta = Sub(scope, Sigmoid(scope, z), labels);
#endif
  debug_show("delta", delta);

  // Calculate xi_0. The README has ξₗᵢ = \sum_{k=0}^{Mₗ-1}(δₗₖ wₖᵢ).
  // Adding an extra dimension for the batch size and leaving l away, that becomes:
  //
  //    ξᵢₛ = \sum_{k=0}^{Mₗ-1}(δₖₛ wₖᵢ) = Wᵀ×δ
  //
  // Not used: we only have a single layer (xi_0 never needs to be calculated).
//  auto xi_0 = MatMul(scope, Transpose(scope, weights_bias, {1, 0}), delta);
//  debug_show("xi_0", xi_0);

  // Update weights_bias. The README has wᵢⱼ' = wᵢⱼ - α δₗᵢ xⱼ, where δₗᵢ xⱼ is the gradient.
  // Adding an extra dimension for the batch size and leaving l away, that becomes:
  //
  //    Gᵢⱼ = \sum_{s=0}^{batch_size-1}(δᵢₛ xⱼₛ) = δ×Xᵀ
  //    wᵢⱼ' = wᵢⱼ - (α / batch_size) Gᵢⱼ
  //
  auto gradient = MatMul(scope.WithOpName("bsG"), delta, Transpose(scope, inputs_1, {1, 0}));
  debug_show("gradient", gradient);
  auto dg = Multiply(scope.WithOpName("alphaG"), learning_rate, gradient);
  debug_show("dg", dg);

  // Update weights_bias.
  auto update_weights_bias_op = AssignSub(scope.WithControlDependencies({update_weights_bias_op2}), weights_bias, dg, AssignSub::Attrs().UseLocking(true)).operation;
  auto updated_weights_bias = Identity(scope.WithControlDependencies({update_weights_bias_op}), weights_bias);

  // Write the graph to a file.
  dump_graph(scope, "neural_network");

  //---------------------------------------------------------------------------
  // Fill training samples.

  // Create a random number generator.
//  std::random_device rd;
//  std::default_random_engine generator(rd());
  std::default_random_engine generator(seed2);
  std::uniform_real_distribution<float> weights_distribution(0.1, 1.0);

  std::vector<float> labels_data; // = { 0, 1, 1, 0, 0, 1, 0, 0, 1 };
  for (int i = 0; i < size_of_batch; ++i)
    labels_data.push_back((weights_distribution(generator) > 0.55) ? 0.f : 1.f);
  int64 const number_of_data_points = labels_data.size();

  // Create a tensor with labels.
  Tensor labels_tensor{DT_FLOAT, TensorShape{number_of_data_points}};
  std::copy_n(labels_data.begin(), labels_data.size(), labels_tensor.flat<float>().data());

  // Create a tensor with inputs.
  Tensor inputs_tensor{DT_FLOAT, TensorShape{2, number_of_data_points}};
  auto tensor_map = inputs_tensor.tensor<float, 2>();
  for (int i = 0; i < number_of_data_points; ++i)
  {
    float x0 = i;
    float rd = weights_distribution(generator);
    float delta_x1 = ((labels_data[i] > 0.5) ? 1 : -1) * rd;    // 0 means below the line (red), 1 means above (green).
    std::cout << "delta_x1 = " << delta_x1 << std::endl;

    tensor_map(0, i) = x0;
    tensor_map(1, i) = x1(x0) + delta_x1;
  }

  std::cout << "labels_tensor is now " << labels_tensor.DebugString(number_of_data_points) << std::endl;
  std::cout << "inputs_tensor is now " << inputs_tensor.DebugString(2 * number_of_data_points) << std::endl;

  // Plot the input data points.
  std::array<std::vector<double>, 2> xp;
  std::array<std::vector<double>, 2> yp;
  for (int i = 0; i < number_of_data_points; ++i)
  {
    int index = std::round(labels_data[i]);
    assert(index == 0 || index == 1);
    xp[index].push_back(tensor_map(0, i));
    yp[index].push_back(tensor_map(1, i));
  }
  assert(xp[0].size() + xp[1].size() == number_of_data_points);

  draw_points(h1, xp, yp);

  //---------------------------------------------------------------------------
  // Running.

  ClientSession session(scope);

  // Fill weights_bias with random initial values.
  // Fill inputs_1 with inputs + ones.
  TF_CHECK_OK(session.Run(ClientSession::FeedType{{inputs, inputs_tensor}}, {}, initialization_operations, nullptr));

  TF_CHECK_OK(session.Run(ClientSession::FeedType{{inputs, inputs_tensor}, {labels, labels_tensor}},
        debug_show.fetch_outputs(), {}, debug_show.outputs()));
  debug_show.dump();
//  debug_show.clear();

  //---------------------------------------------------------------------------
  // Training.

  // Associate remaining placeholder with the input data.
  ClientSession::FeedType inputs_feed = { {labels, labels_tensor} };
  std::vector<Output> fetch_outputs1 = { updated_weights_bias };
  std::vector<Output> fetch_outputs2 = { updated_weights_bias, outputs, loss, residual };

  int plots = 0;
  for (int n = 0; n < 100000; ++n)
  {
    std::vector<Tensor> outputs;
    if (n % 1 != 0)
      TF_CHECK_OK(session.Run(inputs_feed, fetch_outputs1, &outputs));
    else
    {
      std::cout << "seed1 = " << seed1 << "; seed2 = " << seed2 << '\n';
//      debug_show.list();
//      TF_CHECK_OK(session.Run(ClientSession::FeedType{{inputs, inputs_tensor}, {labels, labels_tensor}}, debug_show.fetch_outputs(), {}, debug_show.outputs()));
//      debug_show.dump();
      TF_CHECK_OK(session.Run(inputs_feed, fetch_outputs2, &outputs));
//      Tensor const& wb = debug_show.outputs()->operator[](0);
      auto tensor_map0 = outputs[0].tensor<float, 2>();
      if (!tensor_map0.data())
        continue;
      float xm = tensor_map0(0, 0);
      float ym = tensor_map0(0, 1);
      float om = tensor_map0(0, 2);
      std::cout << "y = " << (-xm / ym) << " * x - " << (om / ym) << std::endl;
      std::cout << "weights_bias = [" << xm << ", " << ym << ", " << om << "]" << std::endl;
      std::ostringstream function_name;
      function_name << "f" << n << "(x)";
      std::ostringstream equation;
      equation << function_name.str() << " = " << (-xm / ym) << " * x - " << (om / ym) << '\n';
#if 1
      auto tensor_map1 = outputs[1].tensor<float, 2>();
      std::cout << "outputs = [";
      char const* sep = "";
      for (int s = 0; s < size_of_batch; ++s)
      {
        std::cout << sep << tensor_map1(0, s);
        sep = ", ";
      }
      std::cout << "]" << std::endl;
      auto loss_tensor_map = outputs[2].tensor<float, 2>();
      std::cout << "loss = [";
      sep = "";
      for (int s = 0; s < size_of_batch; ++s)
      {
        std::cout << sep << std::format("{: >10.7f}", loss_tensor_map(0, s));
        sep = ", ";
      }
      std::cout << "]" << std::endl;
      std::cout << "labels = " << labels_tensor.DebugString(size_of_batch) << std::endl;
      std::cout << "residual = " << outputs[3].DebugString(size_of_batch) << std::endl;
      auto residual_tensor_map = outputs[3].tensor<float, 2>();
      bool correct_prediction = true;
      for (int s = 0; s < size_of_batch; ++s)
      {
        if (std::abs(residual_tensor_map(0, s)) > 0.5f)
        {
          correct_prediction = false;
          break;
        }
      }
      if (correct_prediction)
      {
        std::cout << "Found correct prediction!" << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        break;
      }
#endif
      if (plots++ == 10)
      {
        gnuplot_resetplot(h1);
        plots = 0;
        draw_points(h1, xp, yp);
      }
      gnuplot_cmd(h1, equation.str().c_str());
      gnuplot_plot_equation(h1, function_name.str().c_str(), ("epoch " + std::to_string(n)).c_str());
      std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }
  }

  // Close plot window.
  gnuplot_close(h1);
}
