#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <fstream>
#include <random>

using namespace tensorflow;
using namespace tensorflow::ops;

// Function to create a dense (fully connected) layer.
Output SigmoidDenseLayer(Scope& scope, Input inputs, Input labels, int input_units, int output_units, std::vector<Operation>& need_running,
    std::vector<Output>* debug_outputs)
{
  // Define a weights matrix.
  auto weights = Variable(scope, {output_units, input_units}, DT_FLOAT);

  // Initialization of the weights.
  need_running.push_back(Assign(scope, weights, RandomNormal(scope, {output_units, input_units}, DT_FLOAT)).operation);

  // Define a biases vector.
  auto biases = Variable(scope, {output_units}, DT_FLOAT);

  // Initialization of the biases.
  need_running.push_back(Assign(scope, biases, RandomNormal(scope, {output_units}, DT_FLOAT)).operation);

  // Make a copy of the weights and biases for inspection.
  if (debug_outputs)
  {
    debug_outputs->push_back(Identity(scope, weights));
    debug_outputs->push_back(Identity(scope, biases));
  }

  // Perform `weights * inputs + biases` and apply activation function.
  auto output = Sigmoid(scope, Add(scope, MatMul(scope, weights, inputs), biases));

  // Define a loss function, lets use Mean Squared Error (MSE) (totally random).
  auto loss = ReduceMean(scope, Square(scope, Subtract(scope, labels, output)), {0});

  return output;
}

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

float y(float x)
{
  return 2.5 * x - 4.3;
}

int main()
{
  // Define the scope.
  Scope scope = Scope::NewRootScope();

  //---------------------------------------------------------------------------
  // Create graph.

  // Define the input placeholder with shape [UnknownRank, 2], where 'UnknownRank' means the first dimension is dynamic.
  auto inputs = Placeholder(scope.WithOpName("inputs"), DT_FLOAT, Placeholder::Shape({UnknownRank, 2}));
  // Define the labels placeholder with shape [UnknownRank], where 'UnknownRank' means the dimension is dynamic.
  auto labels = Placeholder(scope.WithOpName("labels"), DT_FLOAT, Placeholder::Shape({UnknownRank}));

  //        .---.
  // x ---> |   | ---> [tanh activation] ---> { >0 --> label = 1.
  // y ---> |   |                             { <0 --> label = 0.
  //        '---'
  //
  // v is the input vector (x, y), or (x1, x2).
  // w is the weights vector (w1, w2).
  // b is the bias, a scalar.
  //
  // Pre-activation output = v ⋅ w + b = x1⋅w1 + x2⋅w2 + b.
  // Since the activation function tanh leaves the sign unchanged, we want the pre-activation output
  // to be zero when a point is exactly on the line, positive when it is above the line, and negative
  // when it is below the line. This means that w2 must be positive.
  //
  // Using the above y(float x): x2 = 2.5 * x1 - 4.3,
  // we have: x1⋅w1 + (2.5 * x1 - 4.3)⋅w2 + b = 0, and w2 > 0 -->
  // x1⋅(w1 + 2.5⋅w2) - 4.3⋅w2 + b = 0 for any x1 -->
  // w1 + 2.5⋅w2 = 0 and -4.3⋅w2 + b = 0 -->
  // w1 = -2.5⋅w2, and b = 4.3⋅w2.
  // For example, if w2=1 > 0, then w1=-2.5 and b=4.3.
  //
  // This proves that a single perceptron with two inputs is able to learn
  // how to separate 2D points from laying above or below a line.

  // Create a neural network layer existing of a single perceptron with two inputs and one output.
  std::vector<Operation> initialization_operations;
  std::vector<Output> debug_outputs;
  auto output = SigmoidDenseLayer(scope, inputs, labels, 2, 1, initialization_operations, &debug_outputs);

  // Gradients of loss w.r.t weights and biases
//  std::vector<tensorflow::Output> grad_outputs;
//  TF_CHECK_OK(AddSymbolicGradients(scope, {loss}, {weights, biases}, &grad_outputs));

  // Assuming you have a learning rate defined as a constant or placeholder.
//  auto optimizer = ApplyGradientDescent(scope, /* trainable variables */, /* learning rate */, /* gradient computation */);

  // Write the graph to a file.
  dump_graph(scope, "neural_network");

  //---------------------------------------------------------------------------
  // Fill training samples.

  // Create a random number generator.
  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.001, 0.1);

  std::vector<float> labels_data = { 0, 1, 1, 0, 0, 1, 0, 0, 1 };
  int64 const number_of_data_points = labels_data.size();
  std::vector<float> inputs_data;
  for (int i = 0; i < number_of_data_points; ++i)
  {
    float x = i;
    float delta_y = (labels_data[i] == 1 ? 1 : -1) * distribution(generator);

    inputs_data.push_back(x);
    inputs_data.push_back(y(x) + delta_y);
  }

  // Create a tensor for inputs and labels.
  Tensor inputs_tensor{DT_FLOAT, TensorShape{number_of_data_points, 2}};
  Tensor labels_tensor{DT_FLOAT, TensorShape{number_of_data_points}};

  // Fill the tensors with our training data.
  std::copy_n(labels_data.begin(), labels_data.size(), labels_tensor.flat<float>().data());
  std::copy_n(inputs_data.begin(), inputs_data.size(), inputs_tensor.flat<float>().data());

  std::cout << "labels_tensor is now " << labels_tensor.DebugString(number_of_data_points) << std::endl;
  std::cout << "inputs_tensor is now " << inputs_tensor.DebugString(number_of_data_points * 2) << std::endl;

  //---------------------------------------------------------------------------
  // Running.

  ClientSession session(scope);

  // Fill weights and biases with random initial values.
  TF_CHECK_OK(session.Run({}, {}, initialization_operations, nullptr));

  // Read back and print debug outputs (initial values of weights and biases).
  std::vector<Tensor> outputs;
  TF_CHECK_OK(session.Run({}, debug_outputs, {}, &outputs));
  std::cout << "Weights: " << outputs[0].DebugString() << std::endl;
  std::cout << "Biases: " << outputs[1].DebugString() << std::endl;

  //---------------------------------------------------------------------------
  // Training.

}
