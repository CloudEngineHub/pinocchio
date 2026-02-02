//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <cmath>
#include <limits>

#include <boost/math/constants/constants.hpp>

#include "pinocchio/macros.hpp"
#include "pinocchio/eigen-common.hpp"
#include "pinocchio/alloca.hpp"

#include "pinocchio/utils/static-if.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/math/comparison-operators.hxx"
#include "pinocchio/math/fwd.hxx"
#include "pinocchio/math/alias.hxx"
#include "pinocchio/math/taylor-expansion.hxx"
#include "pinocchio/math/matrix.hxx"
#include "pinocchio/math/sincos.hxx"
#include "pinocchio/math/rotation.hxx"
#include "pinocchio/math/quaternion.hxx"
#include "pinocchio/math/triangular-matrix.hxx"
#include "pinocchio/math/matrix-block.hxx"
// IWYU pragma: end_exports
