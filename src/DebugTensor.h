#pragma once

#include <array>
#include <vector>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

using FloatType = float;

template<int d_>
class Tensor
{
 public:
  using Index = std::array<int, d_>;    // For example if d=3: i,j,k
  using Shape = Index;                  // Same type :/.
  using Flat = std::vector<FloatType>;

 private:
  Shape shape_;                         // For example, N x M x K where 0 <= i < N, 0 <= j < M, 0 <= k < K.
  Flat flat_;

  int index2fi(Index const& index) const
  {
    int i = 0;
    for (int d = d_ - 1; d >= 0; --d)
    {
      assert(index[d] < shape_[d]);
      i *= shape_[d];
      i += index[d];
    }
    return i;
  }

  int calc_size(int sd) const
  {
    int s = 1;
    for (int d = 0; d < sd; ++d)
      s *= shape_[d];
    return s;
  }

 public:
  template<int d = d_, typename std::enable_if<d == 0, int>::type = 0>
  Tensor(FloatType value = 0) : flat_(1, value) { }

  Tensor(Shape shape, FloatType value = 0) : shape_(shape), flat_(calc_size(d_), value) { assert(!flat_.empty()); }

  FloatType& operator[](Index const& index) { return flat_[index2fi(index)]; }
  FloatType operator[](Index const& index) const { return flat_[index2fi(index)]; }

  template<typename... Args>
  FloatType& operator()(Args&&... args) { return flat_[index2fi({args...})]; }

  template<typename... Args>
  FloatType operator()(Args&&... args) const { return flat_[index2fi({args...})]; }

  Tensor& operator+=(Tensor const& tensor)
  {
    assert(shape_ == tensor.shape_);
    for (int fi = 0; fi < flat_.size(); ++fi)
      flat_[fi] += tensor.flat_[fi];
    return *this;
  }

  friend Tensor operator+(Tensor const& tensor1, Tensor const& tensor2)
  {
    Tensor result(tensor1);
    result += tensor2;
    return result;
  }

  Tensor& operator-=(Tensor const& tensor)
  {
    assert(shape_ == tensor.shape_);
    for (int fi = 0; fi < flat_.size(); ++fi)
      flat_[fi] -= tensor.flat_[fi];
    return *this;
  }

  friend Tensor operator-(Tensor const& tensor1, Tensor const& tensor2)
  {
    Tensor result(tensor1);
    result -= tensor2;
    return result;
  }

  void next_skip(Index& index, int m) const
  {
    for (int d = 0; d < d_; ++d)
    {
      if (d == m)
        continue;
      if (++index[d] < shape_[d])
        return;
      index[d] = 0;
    }
    index[0] = -1;
  }

  Shape const& shape() const { return shape_; }
  int shape(int d) const { return shape_[d]; }

  Flat& flat() { return flat_; }
  Flat const& flat() const { return flat_; }

  // Append one element to the dimension dim.
  // I.e. [AxBxCx...xLx...xK] --> [AxBxCx...x(L+1)x...xK] having appended a [AxBxCx...xK].
  Tensor& append(int dim, Tensor<d_ - 1> const& t)
  {
#ifndef NDEBUG
    // Check that AxBxCx...xK are the same.
    int dt = 0;
    for (int d = 0; d < d_; ++d)
    {
      if (d == dim)
        continue;
      assert(shape_[d] == t.shape(dt++));
    }
#endif
    // Calculate 'AxBxCx...' (all the way up to but not including L)
    int appendsize = 1;
    for (int d = 0; d < dim; ++d)
      appendsize *= shape_[d];
    int blocksize = appendsize * shape_[dim];
    Flat::const_iterator ti = t.flat().begin();
    shape_[dim]++;
    Flat result;
    for (Flat::const_iterator self = flat_.begin(); self != flat_.end(); self += blocksize)
    {
      result.insert(result.end(), self, self + blocksize);
      result.insert(result.end(), ti, ti + appendsize);
      ti += appendsize;
    }
    flat_.swap(result);
    return *this;
  }

  void print_on(std::ostream& os, int fi_offset, int d) const
  {
    if (d == 0) // Scalar
    {
      assert(fi_offset >= 0);
      assert(fi_offset < flat_.size());
      os << flat_[fi_offset];
      return;
    }

    os << '[';
    int offset = fi_offset;
    char const* sep = "";
    int sz = calc_size(d - 1);
    for (int e = 0; e < shape_[d - 1]; ++e)
    {
      os << sep;
      sep = ", ";
      print_on(os, offset, d - 1);
      offset += sz;
    }
    os << ']';
  }

  void print_on(std::ostream& os) const
  {
    print_on(os, 0, d_);
  }

  friend std::ostream& operator<<(std::ostream& os, Tensor const& tensor)
  {
    tensor.print_on(os);
    return os;
  }
};

template<int d_>
bool at_end(std::array<int, d_> const& index)
{
  return index[0] == -1;
}

template<int d1, int d2>
Tensor<d1 + d2 - 2> contract(Tensor<d1> const& tensor1, Tensor<d2> const& tensor2, int m1, int m2)
{
  assert(tensor1.shape(m1) == tensor2.shape(m2));

  std::array<int, d1 + d2 - 2> shape;
  {
    int i = 0;
    for (int i1 = 0; i1 < d1; ++i1)
      if (i1 != m1)
        shape[i++] = tensor1.shape(i1);
    for (int i2 = 0; i2 < d2; ++i2)
      if (i2 != m2)
        shape[i++] = tensor2.shape(i2);
  }
  Tensor<d1 + d2 - 2> result(shape);
  typename Tensor<d1 + d2 - 2>::Index i{};
  for (typename Tensor<d2>::Index i2{}; !at_end<d2>(i2); tensor2.next_skip(i2, m2))
    for (typename Tensor<d1>::Index i1{}; !at_end<d1>(i1); tensor1.next_skip(i1, m1))
    {
      FloatType sum = 0;
//      std::cout << "Setting sum to 0\n";
      for (int dd = 0; dd < tensor1.shape(m1); ++dd)
      {
        i1[m1] = dd;
        i2[m2] = dd;
//        std::cout << "Adding " << tensor1[i1] << " * " << tensor2[i2] << " = " << (tensor1[i1] * tensor2[i2]) << "\n";
        sum += tensor1[i1] * tensor2[i2];
      }
//      std::cout << "Calculated value " << sum << std::endl;
      result[i] = sum;
      result.next_skip(i, -1);
    }
  return result;
}

using Matrix = Tensor<2>;
using Vector = Tensor<1>;
using Scalar = Tensor<0>;

static constexpr FloatType sigmoid(FloatType v)
{
  return FloatType{1} / (FloatType{1} + std::exp(-v));
}

static constexpr FloatType square(FloatType x)
{
  return x * x;
}

class Function
{
 private:
  std::function<FloatType (FloatType)> a_;

 public:
  Function(std::function<FloatType (FloatType)> a) : a_(a) { }

  template<int d_>
  Tensor<d_> operator()(Tensor<d_> const& tensor) const
  {
    Tensor<d_> result(tensor.shape());
    auto const& flat_in = tensor.flat();
    auto& flat_out = result.flat();
    for (int e = 0; e < flat_in.size(); ++e)
      flat_out[e] = a_(flat_in[e]);
    return result;
  }
};

static Vector mse(Matrix const& outputs, Matrix const& targets)
{
  int M = outputs.shape()[0];
  int batch_size = outputs.shape()[1];
  Function sqr(square);
  auto squares = sqr(outputs - targets);
  Vector result{{batch_size}};
  for (int s = 0; s < batch_size; ++s)
  {
    FloatType sum = 0;
    for (int i = 0; i < M; ++i)
      sum += squares(i, s);
    result(s) = sum / M;
  }
  return result;
}

static Matrix derivative_mse(Matrix const& outputs, Matrix const& targets)
{
  int M = outputs.shape()[0];
  return Function([M](FloatType in){ return 2 * in / M; })(outputs - targets);
}
