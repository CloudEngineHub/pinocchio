//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: begin_keep
#include <cstddef>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <vector>
#include <array>

#include <Eigen/Core>

#include "pinocchio/macros.hpp"
#include "pinocchio/eigen-common.hpp"
#include "pinocchio/utils/template-template-parameter.hpp"

#include "pinocchio/math.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports

namespace pinocchio
{
  template<typename T, typename Enable = void>
  struct sizeInBytesImpl;
}

#include "pinocchio/utils/promote-static-eval.hpp"
#include "pinocchio/utils/size-in-bytes.hxx"
#include "pinocchio/utils/eigen.hxx"
#include "pinocchio/utils/std-array.hxx"
#include "pinocchio/utils/std-vector.hxx"
// IWYU pragma: end_exports
