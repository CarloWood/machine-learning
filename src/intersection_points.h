#include <array>
#include <vector>
#include <concepts>
#include <stdexcept>
#include <iostream>

namespace intersections {

// An n-dimensional vector (a column vector with n elements).
template<std::floating_point FloatType, int n>
struct Vector
{
  std::array<FloatType, n> v_;                          // The elements of the vector.

  // Construct a zeroed vector.
  Vector() : v_{} { }

  // Initialize this Vector from an initializer list.
  Vector(std::initializer_list<FloatType> v)
  {
    if (v.size() != n)
      throw std::invalid_argument("Initializer list must have exactly n elements.");
    std::copy(v.begin(), v.end(), v_.begin());
  }

  // Element access.
  FloatType& operator[](int i) { return v_[i]; }
  FloatType operator[](int i) const { return v_[i]; }

  // Add v2.
  Vector& operator+=(Vector const& v2)
  {
    for (int i = 0; i < n; ++i)
      v_[i] += v2[i];
    return *this;
  }

  // Add v1 and v2.
  friend Vector operator+(Vector const& v1, Vector const& v2)
  {
    Vector result(v1);
    result += v2;
    return result;
  }

  // Subtract v2.
  Vector& operator-=(Vector const& v2)
  {
    for (int i = 0; i < n; ++i)
      v_[i] -= v2[i];
    return *this;
  }

  // Subtract v2 from v1.
  friend Vector operator-(Vector const& v1, Vector const& v2)
  {
    Vector result(v1);
    result -= v2;
    return result;
  }

  // Return the dot product of v1 and v2.
  friend FloatType operator*(Vector const& v1, Vector const& v2)
  {
    FloatType result = 0;
    for (int i = 0; i < n; ++i)
      result += v1[i] * v2[i];
    return result;
  }

  // Elementwise multiply with scalar g.
  friend Vector operator*(FloatType g, Vector const& v)
  {
    Vector result(v);
    for (int i = 0; i < n; ++i)
      result[i] *= g;
    return result;
  }

  // For debugging purposes.
  void print_on(std::ostream& os) const
  {
    char const* sep = "";
    os << '(';
    for (int i = 0; i < n; ++i)
    {
      os << sep << v_[i];
      sep = ", ";
    }
    os << ')';
  }

  friend std::ostream& operator<<(std::ostream& os, Vector const& v)
  {
    v.print_on(os);
    return os;
  }
};

// An n-1 dimensional hyper-plane, orthogonal to a given normal vector N, with offset dN from the origin.
template<std::floating_point FloatType, int n>
struct HyperPlane
{
  using VectorType = Vector<FloatType, n>;

  VectorType N_;                                        // The normal of the hyper-plane.
  FloatType d_;                                         // The signed distance (in Normal vectors) from origin to hyper-plane.

  // Create a hyper-plane that satisfies N·X + d = 0.
  HyperPlane(VectorType const& N, FloatType d) : N_(N), d_(d) { }

  // Return the number of N_'s that have to be added to P to end up on this HyperPlane.
  // The sign tells you which "side" of the hyper-plane the point is on.
  FloatType distance(VectorType const& P) const
  {
    FloatType dist = -(N_ * P + d_);
    return dist;
  }

  // Return intersection of the line through C1 and C2 with this HyperPlane.
  VectorType intersection(VectorType const& C1, VectorType const& C2) const
  {
    // Let E be a line through C1 and C2: E: C1 + g(C2 - C1), where g parameterizes the points on E.
    // Fill that in in the line equation to find the intersection:
    // N·(C1 + g(C2 - C1)) + d = 0 --> N·C1 + d + g N·(C2 - C1) = 0 --> g = -(N·C1 + d) / N·(C2 - C1)
    VectorType diff = C2 - C1;
    FloatType g = -(N_ * C1 + d_) / (N_ * diff);
    return C1 + g * diff;
  }
};

template<std::floating_point FloatType, int n>
struct HyperBlock
{
  using VectorType = Vector<FloatType, n>;
  static constexpr int number_of_corners = 1 << n;

  std::array<VectorType, number_of_corners> C_;         // The 2^n corners of the hyper-block.

  // Construct an axis-aligned hyper-block from two adjecent corner vectors.
  HyperBlock(VectorType const& c1, VectorType const& c2)
  {
    VectorType base = c2 - c1;
    for (int ci = 0; ci < number_of_corners; ++ci)
    {
      VectorType c;
      for (int d = 0; d < n; ++d)
      {
        int bit = 1 << d;
        if ((ci & bit))
          c[d] += base[d];
      }
      C_[ci] = c1 + c;
    }
  }

  std::vector<VectorType> intersection_points(HyperPlane<FloatType, n> const& plane);
};

template<std::floating_point FloatType, int n>
std::vector<typename HyperBlock<FloatType, n>::VectorType> HyperBlock<FloatType, n>::intersection_points(HyperPlane<FloatType, n> const& plane)
{
  std::vector<VectorType> intersections;

  std::array<bool, number_of_corners> side;
  for (int ci = 0; ci < number_of_corners; ++ci)
    side[ci] = plane.distance(C_[ci]) <= 0;

  for (int ci = 0; ci < number_of_corners; ++ci)
  {
    for (int d = 0; d < n; ++d)
    {
      int bit = 1 << d;
      int ci2 = ci & ~bit;
      if (side[ci] != side[ci2])
      {
        // Found two corners on opposite sides of the hyper-plane.
        intersections.push_back(plane.intersection(C_[ci], C_[ci2]));
      }
    }
  }

  return intersections;
}

} // namespace intersections
