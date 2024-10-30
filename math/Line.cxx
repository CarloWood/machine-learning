#include "sys.h"
#include "Line.h"
#include "debug.h"

namespace math {

Point Line::intersection_with(Line const& L1) const
{
  // Line0: P0 + λ D0
  // Line1: P1 + ξ D1

  Point const& P0 = point_;
  Direction const& D0 = direction_;
  Point const& P1 = L1.point();
  Direction const& D1 = L1.direction();

  // Let N1 be D1 rotated counter-clockwise by PI/2 (this is floating-point round off error free).
  Direction N1 = D1.normal();

  // Take dot product of D0 with N1:
  double D0_dot_N1 = D0.dot(N1);

  //            intersection
  //                 \ /         P1
  //  --------+-------+-<--------+
  //       ^  |      /    D1  1  |
  //       |  |     /λ           |N1
  //     a |  |_  _/            1|
  //       | ^|   /|             |
  //       |b|| 1/               v
  //       | || /D0
  //       v v|/
  //          +P0
  //         /
  // intersection: P0 + λ D0 = P1 + ξ D1 -->
  //
  //   P1 - P0 = λ D0 - ξ D1
  //
  // Take dot product with N1 (D1 rotated PI/2 radians):
  //
  // λ = (P1 - P0)·N1 / (D0·N1)
  //
  // Note, in the above picture: a = -(P1 - P0)·N1, and b = -D0·N1

  // Take dot product of P1-P0 with N1:
  double P1P0_dot_N1 = (L1.point().x() - point_.x()) * N1.x() + (L1.point().y() - point_.y()) * N1.y();

  // Calculate lambda.
  double lambda = P1P0_dot_N1 / D0_dot_N1;

  // Return intersection point.
  return {point_.x() + lambda * direction_.x(), point_.y() + lambda * direction_.y()};
}

} // namespace math
