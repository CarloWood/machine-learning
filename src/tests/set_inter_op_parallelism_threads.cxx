#include "tensorflow/core/public/session.h"
#include "tensorflow/core/protobuf/config.pb.h"

using namespace tensorflow;

int main()
{
  // Configure session options.
  tensorflow::SessionOptions options;
  options.config.set_inter_op_parallelism_threads(2);

  // Create a new session with the configured options.
  Session* session;
  Status status = NewSession(options, &session);
  if (!status.ok()) {
    std::cerr << status.ToString() << std::endl;
    return 1;
  }

  // Your TensorFlow operations go here.

  // Clean up.
  status = session->Close();
  if (!status.ok()) {
    std::cerr << status.ToString() << std::endl;
    return 1;
  }
}
