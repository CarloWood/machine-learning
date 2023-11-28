#pragma once

#include <string>

namespace tensorflow {
class Scope;
}

void dump_graph(tensorflow::Scope const& scope, std::string name);
