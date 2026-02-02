//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "pinocchio/eigen-common.hpp"

#include "pinocchio/utils/static-if.hpp"

#include "pinocchio/math/sincos.hpp"
#include "pinocchio/math/taylor-expansion.hpp"

#include "pinocchio/context.hxx"
#include "pinocchio/spatial/fwd.hxx"
#include "pinocchio/spatial/motion.hpp"
#include "pinocchio/spatial/skew.hpp"
#include "pinocchio/spatial/se3.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/spatial/log-common.hxx"
#include "pinocchio/spatial/explog.hxx"
#include "pinocchio/spatial/explog-quaternion.hxx"
#include "pinocchio/spatial/log.hxx"
// IWYU pragma: end_exports
