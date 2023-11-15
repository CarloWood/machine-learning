#include <tensorflow/core/platform/env.h>
#include <tensorflow/core/public/session.h>

#include <iostream>

using namespace tensorflow;

int main()
{
  // Configure session options.
  SessionOptions session_options;
  session_options.config.set_inter_op_parallelism_threads(2);

  Session* session;
  Status status = NewSession(session_options, &session);
  if (!status.ok())
  {
    std::cout << status.ToString() << "\n";
    return 1;
  }
  std::cout << "Session successfully created.\n";
}
