#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <iostream>
#include <fstream>

// Basic math operations using tensorflow.

// The sample demonstrates how to
// - create various operations.
// - use placeholder and placeholderWithDefault.
// - use different overloads of Run method.

int main(int argc, char** argv)
{
  using namespace tensorflow;
  using namespace tensorflow::ops;

  // create a root scope
  auto scope = Scope::NewRootScope();

  // We are creating various local scopes so that a new
  // session object is created for all the examples.

  // Configure session options.
  tensorflow::SessionOptions session_options;
  session_options.config.set_inter_op_parallelism_threads(2);

  {
    // An example of doing addition on constants.

    ClientSession session(scope, session_options);

    auto a = Const(scope, 2);
    auto b = Const(scope, 3);

    auto c = Add(scope, a, b);

    std::vector<Tensor> outputs;
    TF_CHECK_OK(session.Run({c}, &outputs));

    // We know that it will be scalar.
    // We can also get the underlying data by calling flat.
    std::cout << "Underlying Scalar value -> " << outputs[0].flat<int>() << std::endl;
  }

  {
    // An example of how to supply a variable (i.e. not a constant)
    // whose value is supplied at the time when we run the session.

    ClientSession session(scope, session_options);

    // We will use Placeholder as the type for our variables.
    auto a = Placeholder(scope, DT_INT32);
    auto b = Placeholder(scope, DT_INT32);

    // Define the add operation that takes the placeholders a and b as inputs.
    auto c = Add(scope, a, b);

    // We now specify the values for our placeholders. Note the way that the
    // Run method is called; it is quite different from the previous example.
    //
    // Here we are using this overload of Run method:
    // Run(FeedType const& inputs, std::vector<Output> const& fetch_outputs,
    //        std::vector<Tensor>* outputs) const;
    //
    // Which takes FeedType (alias of std::unordered_map<Output,
    // Input::Initializer, OutputHash>) as the first argument.
    //
    // Note - In std::unordered_map OutputHash is optional So we just need to
    // supply a map whose key has type "Output" and value that respects
    // Initializer.
    //
    // {a,2} & {b,3} would satisfiy this requirement since type 'a' & 'b'
    // is Output.

    std::vector<Tensor> outputs;
    auto status = session.Run({{{a, 2}, {b, 3}}}, {c}, &outputs);
    TF_CHECK_OK(status);

    // We know that it will be scalar.
    // We can also get the underlying data by calling flat.
    std::cout << "Underlying Scalar value -> " << outputs[0].flat<int>() << std::endl;

    // Create a GraphDef object to hold the graph.
    GraphDef graph_def;
    // Use the scope to serialize the graph.
    Status s = scope.ToGraphDef(&graph_def);

    // Write the GraphDef to a file.
    std::ofstream file("exported_graph.pb", std::ios::out | std::ios::binary);
    if (!graph_def.SerializeToOstream(&file))
      std::cerr << "Failed to write graph to exported_graph.pb." << std::endl;
    else
      std::cout << "Output written to exported_graph.pb; run: `python tb.py` to create logs/ and then `tensorboard --logdir=logs` to view.\n";
  }

  {
    // This is yet another example that makes use of Placeholder however
    // this time we want one of the placeholder to have a default value.
    //
    // In other words, it does not need to be specified during the session
    // execution. if you give a new value it would accept it else would use
    // the default value.

    ClientSession session(scope, session_options);

    // Create an input.
    Input defaultAInput{8};

    // We will use Placeholder as the type for our variables.
    auto a = PlaceholderWithDefault(scope, defaultAInput, PartialTensorShape());
    auto b = Placeholder(scope, DT_INT32);

    // Define the add operation that takes the placeholders a and b as inputs.
    auto c = Add(scope, a, b);

    std::vector<Tensor> outputs;

    // In this Run we are not specifying 'a' so its default value i.e. 8 will be used.
    auto status = session.Run({{{b, 3}}}, {c}, &outputs);
    TF_CHECK_OK(status);

    std::cout
        << "Underlying Scalar value (using default placeholder value [8]) -> "
        << outputs[0].flat<int>() << std::endl;

    // Here we do specify a value for placeholder 'a' i.e. 9.
    status = session.Run({{{a, 9}, {b, 3}}}, {c}, &outputs);
    TF_CHECK_OK(status);

    std::cout << "Underlying Scalar value (after supplying new value [9]) -> "
              << outputs[0].flat<int>() << std::endl;
  }
}
