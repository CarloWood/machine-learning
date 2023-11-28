#include "dump_graph.h"
#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <iostream>
#include <fstream>

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
