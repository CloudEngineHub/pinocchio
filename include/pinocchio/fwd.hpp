//
// Copyright (c) 2018-2018 CNRS
// Copyright (c) 2018-2026 INRIA
//

#ifndef __pinocchio_fwd_hpp__
#define __pinocchio_fwd_hpp__

// Forward declaration of the main pinocchio namespace
namespace pinocchio
{
}

#ifdef _WIN32
  #include <windows.h>
  #undef far
  #undef near
#endif

#include <cassert>

#ifdef PINOCCHIO_EIGEN_CHECK_MALLOC
  #ifndef EIGEN_RUNTIME_NO_MALLOC
    #define EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
    #define EIGEN_RUNTIME_NO_MALLOC
  #endif
#endif

#include "pinocchio/macros.hpp"
#include "pinocchio/deprecated.hpp"
#include "pinocchio/warning.hpp"
#include "pinocchio/config.hpp"
#include "pinocchio/unsupported.hpp"

// Include Eigen components
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
  #include <Eigen/AccelerateSupport>
#endif

#include "pinocchio/eigen-macros.hpp"
#include <unsupported/Eigen/CXX11/Tensor>

#include "pinocchio/utils/helpers.hpp"
#include "pinocchio/utils/cast.hpp"
#include "pinocchio/utils/check.hpp"

#include "pinocchio/container/boost-container-limits.hpp"

#ifdef PINOCCHIO_EIGEN_CHECK_MALLOC
  #ifndef EIGEN_RUNTIME_NO_MALLOC
    #define EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
    #define EIGEN_RUNTIME_NO_MALLOC
  #endif
#endif

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
  #include <Eigen/AccelerateSupport>
#endif

#ifdef PINOCCHIO_EIGEN_CHECK_MALLOC
  #ifdef EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
    #undef EIGEN_RUNTIME_NO_MALLOC
    #undef EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
  #endif
#endif

#ifdef PINOCCHIO_WITH_HPP_FCL
  #pragma message PINOCCHIO_WARN(                                                                  \
    "PINOCCHIO_WITH_HPP_FCL define is deprecated, please use PINOCHIO_WITH_COLLISION instead")
  #define PINOCCHIO_WITH_COLLISION
#endif

#include "pinocchio/core/binary-op.hpp"
#include "pinocchio/core/unary-op.hpp"

#include <cstddef> // std::size_t

#include <type_traits>
#include "pinocchio/traits.hxx" // IWYU pragma: export

#include "pinocchio/alloca.hpp"
#include "pinocchio/context.hpp"

#endif // #ifndef __pinocchio_fwd_hpp__
