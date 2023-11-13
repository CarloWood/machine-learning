#include <tensorflow/core/platform/env.h>
#include <tensorflow/core/public/session.h>

#include <iostream>

using namespace tensorflow;

int main()
{
  Session* session;
  Status status = NewSession(SessionOptions(), &session);
  if (!status.ok())
  {
    std::cout << status.ToString() << "\n";
    return 1;
  }
  std::cout << "Session successfully created.\n";
}
