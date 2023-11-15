#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>

__asm__(".symver omp_in_parallel,omp_in_parallel@OMP_1.0");

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
  TF_CHECK_OK(session.Run({joinOp}, &outputs));

  // See our output using DebugString that tells more information about the tensor.
  std::cout << "DebugString -> " << outputs[0].DebugString() << std::endl;

  // We can also get the underlying data by calling flat.
  std::cout << "Underlying Scalar value -> " << outputs[0].flat<tstring>() << std::endl;
}
