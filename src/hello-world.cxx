#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>

// Simple hello world using TensorFlow.

// The sample demonstrates how to
// - create various ops (Const & StringJoin).
// - pass them around for e.g. StringJoin takes an list of other ops.
// - pass the final ops to the session.
// - get the result of the session.
// - a simple peek inside the output using the DebugString & by flattening it.

int main(int argc, char **argv)
{
  using namespace tensorflow;
  using namespace tensorflow::ops;

  // Create a root scope.
  auto scope = Scope::NewRootScope();

  // Define various constans/inputs on which we will perform an operation.
  auto hello = Const(scope, std::string("hello"));
  auto space = Const(scope, std::string(" "));
  auto world = Const(scope, std::string("world !"));

  // StringJoin operation.
  auto joinOp = StringJoin(scope, {hello, space, world});

  // Configure session options.
  tensorflow::SessionOptions session_options;
  session_options.config.set_inter_op_parallelism_threads(2);

  // Create a session that takes our scope as the root scope.
  ClientSession session(scope, session_options);

  // Run.
  std::vector<Tensor> outputs;
  auto status = session.Run({joinOp}, &outputs);
  TF_CHECK_OK(status);

  // See our output using DebugString that tells more information about the tensor.
  for (Tensor const& tensor : outputs)
    std::cout << "DebugString -> " << tensor.DebugString() << std::endl;

  // We can also get the underlying data by calling flat.
  for (Tensor const& tensor : outputs)
    std::cout << "Underlying Scalar value -> " << tensor.flat<tstring>() << std::endl;
}
