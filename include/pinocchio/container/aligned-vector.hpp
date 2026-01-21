//
// Copyright (c) 2016-2018 CNRS
// Copyright (c) 2018-2025 INRIA
//

#ifndef __pinocchio_container_aligned_vector_hpp__
#define __pinocchio_container_aligned_vector_hpp__

#include "pinocchio/deprecated.hpp"
#include "pinocchio/warning.hpp"

#include <vector>

#pragma message PINOCCHIO_WARN(                                                                    \
  "<pinocchio/container/aligned-vector.hpp> is deprecated and will be removed")

// This macro is deprecated
#define PINOCCHIO_ALIGNED_STD_VECTOR(Type) ::pinocchio::container::aligned_vector<Type>

#define PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(T) PINOCCHIO_ALIGNED_STD_VECTOR(T)

namespace pinocchio
{
  namespace container
  {

    template<typename T>
    using aligned_vector PINOCCHIO_DEPRECATED_MESSAGE(
      "aligned_vector is deprecated, please use std::vector") = std::vector<T>;

  } // namespace container

} // namespace pinocchio

#endif // ifndef __pinocchio_container_aligned_vector_hpp__
