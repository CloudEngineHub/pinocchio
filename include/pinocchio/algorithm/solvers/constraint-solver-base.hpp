//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: begin_keep
#include <cstddef>
#include <limits>
#include <vector>

#include "pinocchio/macros.hpp"
#include "pinocchio/fwd.hpp"
#include "pinocchio/utils/check.hpp"

#include "pinocchio/math.hpp"

#ifdef PINOCCHIO_WITH_COLLISION
  #include <coal/timings.h>
#endif // PINOCCHIO_WITH_COLLISION
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/algorithm/solvers/constraint-solver-base.hxx"
// IWYU pragma: end_exports