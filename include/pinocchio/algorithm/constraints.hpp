//
// Copyright (c) INRIA 2026
//

// IWYU pragma: begin_keep
#include <cassert>

#include <boost/blank.hpp>

#include <Eigen/Core>

#include "pinocchio/fwd.hpp"
#include "pinocchio/macros.hpp"
#include "pinocchio/eigen-common.hpp"
#include "pinocchio/context.hpp"

#include "pinocchio/utils/check.hpp"
#include "pinocchio/utils/static-if.hpp"

#include "pinocchio/math.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/algorithm/constraints/fwd.hxx"

#include "pinocchio/algorithm/constraints/sets/set-base.hxx"
#include "pinocchio/algorithm/constraints/sets/box-set.hxx"

#include "pinocchio/algorithm/constraints/sets/cone-base.hxx"
#include "pinocchio/algorithm/constraints/sets/coulomb-friction-cone.hxx"
#include "pinocchio/algorithm/constraints/sets/zero-cone.hxx"
#include "pinocchio/algorithm/constraints/sets/zero-cone-jordan-operation.hxx"
#include "pinocchio/algorithm/constraints/sets/second-order-cone-jordan-operation.hxx"
#include "pinocchio/algorithm/constraints/sets/orthant-cone.hxx"
#include "pinocchio/algorithm/constraints/sets/orthant-cone-jordan-operation.hxx"
#include "pinocchio/algorithm/constraints/sets/full-space-cone.hxx"
// IWYU pragma: end_exports
