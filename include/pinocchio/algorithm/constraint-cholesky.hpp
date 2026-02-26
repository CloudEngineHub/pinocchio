//
// Copyright (c) INRIA 2026
//

#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <vector>

#include <Eigen/Core>

#include "pinocchio/macros.hpp"
#include "pinocchio/deprecated.hpp"
#include "pinocchio/eigen-common.hpp"

#include "pinocchio/utils/alloca.hpp"
#include "pinocchio/utils/check.hpp"
#include "pinocchio/utils/reference.hpp"
#include "pinocchio/utils/size-in-bytes.hpp"

#include "pinocchio/container/eigen-storage.hpp"

#include "pinocchio/math.hpp"
#include "pinocchio/multibody.hpp"
#include "pinocchio/algorithm/constraints.hpp"
#include "pinocchio/algorithm/check-model.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraint-cholesky-decl.hxx"
#include "pinocchio/algorithm/constraint-cholesky-def.hxx"
// IWYU pragma: end_exports
