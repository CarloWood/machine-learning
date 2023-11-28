#include "DebugShow.h"
#include "dump_graph.h"
#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <iostream>

using namespace tensorflow;

int main()
{
  // Define the scope.
  Scope scope = Scope::NewRootScope().ExitOnError();
  DebugShow debug_show(scope);
  std::vector<Operation> initialization_operations;

  //---------------------------------------------------------------------------
  // Graph construction.

  auto x = ops::Placeholder(scope.WithOpName("x"), DT_FLOAT, ops::Placeholder::Shape({2}));
  debug_show("x", x);

  //---------------------------------------------------------------------------
  // Create graph.

  // Get rank 1 vector of size 1 with the element x_1. It is not possible to get this as a scalar!
  auto x_1 = ops::Slice(scope.WithOpName("x_1"), x, {1}, {1});
  // Turn it into a rank 1 vector of size 1 that is true if x_1 < 0.
  auto pred_rank1 = ops::Less(scope.WithOpName("pred"), x_1, {0.f});
  // Now reduce the rank to 0 (a scalar) because ops::Switch demands a rank 0 tensor.
  auto pred = ops::ReduceAny(scope, pred_rank1, 0);

  // If x_1 is negative, multiply x with -1 and that into negative_x, otherwise move it into positive_x.
  auto switch_op = ops::Switch(scope, x, pred);
  Output positive_x = switch_op.output_false;
  Output negative_x = ops::Multiply(scope, -1.f, switch_op.output_true);

  // Which ever one it was moved into now has a positive x[1] element. Move that into abs_x.
  auto abs_x = ops::Merge(scope.WithOpName("abs_x"), { positive_x, negative_x });
  debug_show("abs_x", abs_x.output);

  dump_graph(scope, "switchmerge");

  // Prepare feed data.
  Tensor x_tensor{DT_FLOAT, TensorShape{2}};
  float* x_data = x_tensor.flat<float>().data();
  x_data[0] = 1.0;
  x_data[1] = -2.1;     // x_1 is negative, therefore abs_x will become -x.

  std::cout << "x_tensor = " << x_tensor << '\n';

  //---------------------------------------------------------------------------
  // Running.

  ClientSession session(scope);

  ClientSession::FeedType inputs_feed = { { x, x_tensor } };

  // Run any initialization operations.
  if (!initialization_operations.empty())
    TF_CHECK_OK(session.Run(inputs_feed, {}, initialization_operations, nullptr));

  TF_CHECK_OK(session.Run(inputs_feed, debug_show.fetch_outputs(), {}, debug_show.outputs()));
  debug_show.dump();
}
