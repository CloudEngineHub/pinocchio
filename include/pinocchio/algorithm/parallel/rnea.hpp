//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: begin_keep
#include <cstddef>

#include <omp.h>

#include <Eigen/Core>

#include "pinocchio/macros.hpp"

#include "pinocchio/utils/openmp.hpp"
#include "pinocchio/multibody/pool.hpp"
#include "pinocchio/algorithm/rnea.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/algorithm/parallel/rnea.hxx"
// IWYU pragma: end_exports
