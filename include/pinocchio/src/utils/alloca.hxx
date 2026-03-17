//
// Copyright (c) 2024-2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/utils/alloca.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/utils/alloca.hpp"
#endif // PINOCCHIO_LSP

#define PINOCCHIO_ALLOCA EIGEN_ALLOCA
#define PINOCCHIO_ALIGNED_PTR(ptr, align)                                                          \
  reinterpret_cast<void *>(((intptr_t)ptr + (align - 1)) & ~(align - 1))
#define PINOCCHIO_EIGEN_MAP_ALLOCA(S, rows, cols)                                                  \
  PINOCCHIO_EIGEN_MAP_ALLOCA_ALIGNED(S, rows, cols, EIGEN_DEFAULT_ALIGN_BYTES)
#define PINOCCHIO_EIGEN_MAP_ALLOCA_ALIGNED(S, rows, cols, align)                                   \
  static_cast<S *>(PINOCCHIO_ALIGNED_PTR(                                                          \
    PINOCCHIO_ALLOCA(size_t(rows * cols) * sizeof(S) + (align > 0 ? (align - 1) : 0)), align)),    \
    rows, cols
