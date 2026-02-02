//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <limits>
#include <algorithm>
#include <string>
#include <ostream>
#include <iterator>
#include <sstream>
#include <utility>
#include <vector>
#include <map>
#include <set>

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include <boost/next_prior.hpp>

#include "pinocchio/context.hxx"
#include "pinocchio/traits.hpp"
#include "pinocchio/macros.hpp"
#include "pinocchio/deprecated.hpp"

#include "pinocchio/utils/static-if.hpp"

#include "pinocchio/container/double-entry-container.hpp"
#include "pinocchio/container/matrix-stack.hpp"

#include "pinocchio/common.hpp"
#include "pinocchio/math.hpp"
#include "pinocchio/spatial.hpp"
#include "pinocchio/multibody/joint.hpp"
#include "pinocchio/multibody/liegroup.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/multibody/fwd.hxx"
#include "pinocchio/multibody/model-item.hxx"
#include "pinocchio/multibody/frame.hxx"
#include "pinocchio/multibody/model.hxx"
// IWYU pragma: end_exports
