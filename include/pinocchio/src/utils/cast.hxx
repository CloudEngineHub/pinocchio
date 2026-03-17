//
// Copyright (c) 2020 INRIA
//

#ifndef __pinocchio_utils_cast_hxx__
#define __pinocchio_utils_cast_hxx__

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include <Eigen/Core>
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<typename NewScalar, typename Scalar>
  NewScalar cast(const Scalar & value)
  {
    return Eigen::internal::cast_impl<Scalar, NewScalar>::run(value);
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_utils_cast_hxx__
