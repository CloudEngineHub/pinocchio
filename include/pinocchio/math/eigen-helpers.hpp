//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_eigen_helpers_hpp__
#define __pinocchio_math_eigen_helpers_hpp__

#include "pinocchio/math/fwd.hpp"

namespace pinocchio
{

#define PINOCCHIO_EIGEN_HELPER_FUNC(method)                                                        \
  template<typename Matrix>                                                                        \
  void method(const Eigen::MatrixBase<Matrix> & mat)                                               \
  {                                                                                                \
    mat.const_cast_derived().method();                                                             \
  }

  PINOCCHIO_EIGEN_HELPER_FUNC(setZero);
  PINOCCHIO_EIGEN_HELPER_FUNC(setOnes);
  PINOCCHIO_EIGEN_HELPER_FUNC(setIdentity);

#undef PINOCCHIO_EIGEN_HELPER_FUNC

  constexpr Eigen::AlignmentType to_eigen_alignment(std::size_t alignment_value)
  {
    using Eigen::Aligned;
    using Eigen::Aligned128;
    using Eigen::Aligned16;
    using Eigen::Aligned32;
    using Eigen::Aligned64;
    using Eigen::Aligned8;
    using Eigen::Unaligned;

    switch (alignment_value)
    {
    case 8:
      return Aligned8;
    case 16:
      return Aligned16;
    case 32:
      return Aligned32;
    case 64:
      return Aligned64;
    case 128:
      return Aligned128;
    case 0:
      return Unaligned; // treat 0 as simply unaligned
    default:
      return Aligned; // "default" Eigen alignment policy
    }
  }
} // namespace pinocchio

} // namespace pinocchio

#endif // ifndef __pinocchio_math_eigen_helpers_hpp__
