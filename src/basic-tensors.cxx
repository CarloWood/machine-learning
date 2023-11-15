#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <iostream>

// Basic tensor operations using tensorflow.

// The sample demonstrates how to
// - create various operations.
// - use placeholder and placeholderWithDefault.
// - use different overloads of Run method.

int main(int argc, char** argv)
{
  using namespace tensorflow;
  using namespace tensorflow::ops;

  // Create a root scope.
  auto scope = Scope::NewRootScope();

  {
    // We already know how to create Scalar which is a tensor of Rank 0.
    auto aScalar = Input(2);
    std::cout << "Dimensions of a scalar - " << aScalar.tensor().shape().dims() << std::endl;

    // A tensor of Rank 1 is called a vector.
    auto aVector = Input({2, 3});
    std::cout << "Dimensions of a vector - " << aVector.tensor().shape().dims() << std::endl;

    // A tensor of Rank 2 is called a matrix.
    auto aMatrix = Input({{2, 3}, {6, 5}});
    std::cout << "Dimensions of a matrix - " << aMatrix.tensor().shape().dims() << std::endl;

    // A tensor of Rank 3 or more is not known by any special name tensor.
  }

  // Configure session options.
  tensorflow::SessionOptions session_options;
  // To suppress the Info message "I tensorflow/core/common_runtime/process_util.cc:146]
  // Creating new thread pool with default inter op setting: 2. Tune using inter_op_parallelism_threads for best performance."
  session_options.config.set_inter_op_parallelism_threads(2);
  // Necessary in order to change the graph after running it and then run it again.
  session_options.config.mutable_experimental()->set_disable_optimize_for_static_graph(true);

  {
    // When building the ops you can specify the shape explicitly as well.

    // 2x2 matrix with all elements = 10.
    auto c1 = Const(scope, 10, /* shape */ {2, 2});

    // [1 1] * [41; 1]
    auto x = MatMul(scope, {{1, 1}}, {{41}, {1}});

    ClientSession session(scope, session_options);

    std::vector<Tensor> outputs;

    auto status = session.Run({x}, &outputs);
    TF_CHECK_OK(status);

    std::cout << "Underlying Scalar value -> " << outputs[0].flat<int>() << std::endl;

    // [1 2 3 4] + 10
    auto y = Add(scope, {1, 2, 3, 4}, 10);

    status = session.Run({y}, &outputs);        // <-- Returns an error condition.
    TF_CHECK_OK(status);                        // <-- Error is printed here.

    std::cout << "Underlying vector value -> " << outputs[0].flat<int>() << std::endl;
  }
}
