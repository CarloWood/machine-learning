#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <fstream>
#include <random>

using namespace tensorflow;
using namespace tensorflow::ops;

void dump_graph(Scope const& scope, std::string name)
{
  // Create a GraphDef object to hold the graph.
  GraphDef graph_def;
  // Use the scope to serialize the graph.
  Status s = scope.ToGraphDef(&graph_def);

  // Write the GraphDef to a file.
  std::string filename = name + ".pb";
  std::ofstream file(filename , std::ios::out | std::ios::binary);
  if (!graph_def.SerializeToOstream(&file))
    std::cerr << "Failed to write graph to " << filename << "." << std::endl;
  else
    std::cout << "Output written to \"" << filename << "\"; "
      "run: `python tb.py " << filename << "` to create logs/ and then `tensorboard --logdir=logs` to view.\n";
}

float x1(float x0)
{
  return 2.5f * x0 - 4.3f;
}

class DebugShow
{
 private:
  Scope& scope_;
  std::vector<std::string> labels_;
  std::vector<Output> fetch_outputs_;
  std::vector<Tensor> outputs_;

 public:
  DebugShow(Scope& scope) : scope_(scope) { }

  std::vector<Output> const& fetch_outputs() const { return fetch_outputs_; }
  std::vector<Tensor>* outputs() { return &outputs_; }

  void operator()(std::string label, Output const& output)
  {
    labels_.push_back(label);
    fetch_outputs_.push_back(Identity(scope_, output));
  }

  void operator()(Scope const& scope, std::string label, Output const& output)
  {
    labels_.push_back(label);
    fetch_outputs_.push_back(Identity(scope.WithOpName(label.c_str()), output));
  }

  void dump()
  {
    for (int i = 0; i < labels_.size(); ++i)
      std::cout << labels_[i] << ": " << outputs_[i].DebugString(40) << std::endl;
    outputs_.clear();
  }

  void remove(std::string label)
  {
    for (int i = 0; i < labels_.size(); ++i)
    {
      if (labels_[i] == label)
      {
        labels_.erase(labels_.begin() + i);
        fetch_outputs_.erase(fetch_outputs_.begin() + i);
        return;
      }
    }
  }
};

int main()
{
  // Define the scope.
  Scope scope = Scope::NewRootScope();
  DebugShow debug_show(scope);
  std::vector<Operation> run_operations;

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
  constexpr float alpha = 1;

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

  //---------------------
  // Initialization

  // Initialization of the weights and bias.
  run_operations.push_back(Assign(scope.WithOpName("assignW"), weights_bias,
        RandomNormal(scope.WithOpName("Rand"), {output_units, input_units + 1}, DT_FLOAT, RandomNormal::Attrs().Seed(1))).operation);
  debug_show("weights_bias", weights_bias);

  // Append a 1 to all inputs.
  auto batch_size = Slice(scope.WithOpName("BatchSize"), Shape(scope, inputs), {1}, {1});
  auto inputs_1 = Concat(scope.WithOpName("Inputs1"), {Input(inputs), ExpandDims(scope, Fill(scope, batch_size, 1.0f), 0)}, 0);
  debug_show("inputs_1", inputs_1);

  // Divide alpha by batch_size for use during backpropagation.
  auto learning_rate = Div(scope.WithOpName("Rate"), alpha, Cast(scope, batch_size, DT_FLOAT));
  debug_show("learning_rate", learning_rate);

  //---------------------
  // Forward propagation.

  // Perform `weights_bias * inputs_1` and apply activation function.
  auto outputs = Sigmoid(scope, MatMul(scope, weights_bias, inputs_1));
  debug_show("outputs", outputs);

  // Calculate the difference between the outputs and the targets.
  auto residual = Subtract(scope.WithOpName("Residual"), outputs, labels);     // Shape: [output_units x batch_size]

  // The loss function is not really used.
  {
    // Define a loss function, lets use Mean Squared Error (MSE) (totally random).
    auto loss = ReduceMean(scope.WithOpName("Loss"), Square(scope, residual), {0});
    debug_show("loss", loss);
  }

  //---------------------
  // Backward propagation.

  // Calculate the initial xi value. This is ∂L/∂oᵢ, where L is the loss and O the output (weights_bias * inputs_1).
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
  //
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
  auto update_weights_bias = AssignSub(scope, weights_bias, dg, AssignSub::Attrs().UseLocking(true)).operation;
  debug_show(scope.WithControlDependencies({update_weights_bias}), "weights_bias", weights_bias);

  // Write the graph to a file.
  dump_graph(scope, "neural_network");

  //---------------------------------------------------------------------------
  // Fill training samples.

  // Create a random number generator.
  std::default_random_engine generator;
  generator.seed(1);
  std::uniform_real_distribution<float> distribution(0.001, 0.1);

  std::vector<float> labels_data = { 0, 1, 1, 0, 0, 1, 0, 0, 1 };
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
    float delta_x1 = (labels_data[i] == 1 ? 1 : -1) * distribution(generator);

    tensor_map(0, i) = x0;
    tensor_map(1, i) = x1(x0) + delta_x1;
  }

  std::cout << "labels_tensor is now " << labels_tensor.DebugString(number_of_data_points) << std::endl;
  std::cout << "inputs_tensor is now " << inputs_tensor.DebugString(2 * number_of_data_points) << std::endl;

  //---------------------------------------------------------------------------
  // Running.

  ClientSession session(scope);

  // Fill weights and biases with random initial values.
  TF_CHECK_OK(session.Run({}, {}, run_operations, nullptr));

  // Associate the placeholder with the input data.
  ClientSession::FeedType inputs_feed = {
    {inputs, inputs_tensor},
    {labels, labels_tensor}
  };

  // Run the session with the feed and debug fetches.
  TF_CHECK_OK(session.Run(inputs_feed, debug_show.fetch_outputs(), debug_show.outputs()));
  debug_show.dump();

  //---------------------------------------------------------------------------
  // Training.

  debug_show.remove("inputs");
  debug_show.remove("labels");
  debug_show.remove("inputs_1");
  debug_show.remove("learning_rate");
  debug_show.remove("xi_1");
  debug_show.remove("delta");
  debug_show.remove("gradient");

  debug_show.remove("loss");
  debug_show.remove("weights_bias");
  debug_show.remove("dg");

  for (int n = 0; n < 1000000; ++n)
  {
    // Run the session with the feed and debug fetches.
    TF_CHECK_OK(session.Run(inputs_feed, debug_show.fetch_outputs(), debug_show.outputs()));
    debug_show.dump();
  }
}
