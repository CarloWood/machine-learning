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

constexpr int UnknownRank = -1;

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

  void dump() const
  {
    for (int i = 0; i < labels_.size(); ++i)
      std::cout << labels_[i] << ": " << outputs_[i].DebugString(40) << std::endl;
  }
};

int main()
{
  // Define the scope.
  Scope scope = Scope::NewRootScope();
  DebugShow debug_show(scope);

  //---------------------------------------------------------------------------
  // Create graph.

  // Define the input placeholder with shape [2, UnknownRank], where 'UnknownRank' means the first dimension is dynamic.
  auto inputs = Placeholder(scope.WithOpName("inputs"), DT_FLOAT, Placeholder::Shape({2, UnknownRank}));
  debug_show("inputs", inputs);

  // Define the labels placeholder with shape [UnknownRank], where 'UnknownRank' means the dimension is dynamic.
  auto labels = Placeholder(scope.WithOpName("labels"), DT_FLOAT, Placeholder::Shape({UnknownRank}));
  debug_show("labels", labels);

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
  std::vector<Operation> initialization_operations;
  int input_units = 2;
  int output_units = 1;

  // Define a weights_bias matrix.
  auto weights_bias = Variable(scope, {output_units, input_units + 1}, DT_FLOAT);

  // Initialization of the weights and bias.
  initialization_operations.push_back(Assign(scope, weights_bias, RandomNormal(scope, {output_units, input_units + 1}, DT_FLOAT)).operation);
  debug_show("weights_bias", weights_bias);

  // Get the dynamic shape of 'inputs'.
  auto inputs_shape = Shape(scope, inputs);

  // Extract the batch size (the second dimension of 'inputs_shape')
  auto batch_size = Slice(scope, inputs_shape, {1}, {1});

  // Create a tensor of ones with dynamic shape [batch_size].
  auto ones = Fill(scope, batch_size, 1.0f);
  debug_show("ones", ones);

  // For Concat the work, 'ones' needs to have the shape [batch_size, 1].
  auto ones_expanded = ExpandDims(scope, ones, 0);
  debug_show("ones_expanded", ones_expanded);

  // Concatenate 'inputs' and 'ones_expanded' along the second dimension.
  auto inputs_1 = Concat(scope, {Input(inputs), ones_expanded}, 0);
  debug_show("inputs_1", inputs_1);

  // Perform `weights_bias * inputs_1` and apply activation function.
  auto output = Sigmoid(scope, MatMul(scope, weights_bias, inputs_1));
  debug_show("output", output);

  // Define a loss function, lets use Mean Squared Error (MSE) (totally random).
  auto loss = ReduceMean(scope, Square(scope, Subtract(scope, output, labels)), {0});

  // Write the graph to a file.
  dump_graph(scope, "neural_network");

  //---------------------------------------------------------------------------
  // Fill training samples.

  // Create a random number generator.
  std::default_random_engine generator;
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
  TF_CHECK_OK(session.Run({}, {}, initialization_operations, nullptr));

  // Associate the placeholder with the input data.
  ClientSession::FeedType inputs_feed = {
    {inputs, inputs_tensor},
    {labels, labels_tensor}
  };
  TF_CHECK_OK(session.Run(inputs_feed, debug_show.fetch_outputs(), debug_show.outputs()));
  debug_show.dump();

  //---------------------------------------------------------------------------
  // Training.
}
