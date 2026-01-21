//
// Copyright (c) 2019-2020 INRIA
//

#ifndef __pinocchio_math_tensor_hpp__
#define __pinocchio_math_tensor_hpp__

#include "pinocchio/fwd.hpp"
// #incleade

namespace pinocchio
{

  // Use the default Eigen::Tensor module
  template<typename Scalar_, int NumIndices_, int Options_ = 0, typename IndexType = Eigen::Index>
  using Tensor = Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>;

} // namespace pinocchio

#endif // ifndef __pinocchio_math_tensor_hpp__
