#pragma once

#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <vector>
#include <string>
#include <iostream>

class DebugShow
{
 private:
  tensorflow::Scope& scope_;
  std::vector<std::string> labels_;
  std::vector<tensorflow::Output> fetch_outputs_;
  std::vector<std::string> inactive_labels_;
  std::vector<tensorflow::Output> inactive_fetch_outputs_;
  std::vector<tensorflow::Tensor> outputs_;

 public:
  DebugShow(tensorflow::Scope& scope) : scope_(scope) { }

  std::vector<tensorflow::Output> const& fetch_outputs() const { return fetch_outputs_; }
  std::vector<tensorflow::Tensor>* outputs() { return &outputs_; }

  void operator()(std::string label, tensorflow::Output const& output)
  {
    labels_.push_back(label);
    fetch_outputs_.push_back(tensorflow::ops::Identity(scope_, output));
  }

  void operator()(tensorflow::Scope const& scope, std::string label, tensorflow::Output const& output)
  {
    labels_.push_back(label);
    fetch_outputs_.push_back(tensorflow::ops::Identity(scope.WithOpName(label.c_str()), output));
  }

  void list()
  {
    for (int i = 0; i < labels_.size(); ++i)
      std::cout << i << " : " << labels_[i] << std::endl;
  }

  void dump()
  {
    for (int i = 0; i < labels_.size(); ++i)
      std::cout << labels_[i] << ": " << outputs_[i].DebugString(40) << std::endl;
    outputs_.clear();
  }

  void remove(std::string label)
  {
    for (int i = 0; i < labels_.size(); ++i)
    {
      if (labels_[i] == label)
      {
        inactive_labels_.push_back(label);
        inactive_fetch_outputs_.push_back(fetch_outputs_[i]);
        labels_.erase(labels_.begin() + i);
        fetch_outputs_.erase(fetch_outputs_.begin() + i);
        return;
      }
    }
  }

  void add(std::string label)
  {
    for (int i = 0; i < inactive_labels_.size(); ++i)
    {
      if (inactive_labels_[i] == label)
      {
        labels_.push_back(label);
        fetch_outputs_.push_back(inactive_fetch_outputs_[i]);
        inactive_labels_.erase(inactive_labels_.begin() + i);
        inactive_fetch_outputs_.erase(inactive_fetch_outputs_.begin() + i);
        return;
      }
    }
  }

  void clear()
  {
    for (int i = 0; i < labels_.size(); ++i)
    {
      inactive_labels_.push_back(labels_[i]);
      inactive_fetch_outputs_.push_back(fetch_outputs_[i]);
    }
    labels_.clear();
    fetch_outputs_.clear();
  }
};

