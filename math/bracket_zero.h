#pragma once

#include <functional>

namespace math {

double bracket_zero(double left, double right, std::function<double(double)> const& f);

} // namespace math
