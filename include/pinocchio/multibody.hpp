//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <string>
#include <ostream>
#include <iterator>
#include <sstream>
#include <utility>
#include <memory>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>

#include <boost/next_prior.hpp>
#include <boost/variant/variant.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "pinocchio/context.hpp"
#include "pinocchio/fwd.hpp"
#include "pinocchio/macros.hpp"
#include "pinocchio/deprecated.hpp"
#include "pinocchio/config.hpp"

#include "pinocchio/utils/static-if.hpp"
#include "pinocchio/utils/shared-ptr.hpp"

#include "pinocchio/container/double-entry-container.hpp"
#include "pinocchio/container/matrix-stack.hpp"

#include "pinocchio/common.hpp"
#include "pinocchio/math.hpp"
#include "pinocchio/spatial.hpp"
#include "pinocchio/multibody/joint.hpp"
#include "pinocchio/multibody/liegroup.hpp"
#include "pinocchio/algorithm/fwd.hpp"

#include "pinocchio/serialization/serializable.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/multibody/fwd.hxx"
#include "pinocchio/multibody/force-set.hxx"
#include "pinocchio/multibody/model-item.hxx"
#include "pinocchio/multibody/frame.hxx"
#include "pinocchio/multibody/model.hxx"
#include "pinocchio/multibody/data.hxx"

#include "pinocchio/multibody/coal.hxx"
#include "pinocchio/multibody/instance-filter.hxx"
#include "pinocchio/multibody/geometry-object.hxx"
#include "pinocchio/multibody/geometry.hxx"
#include "pinocchio/multibody/geometry-object-filter.hxx"
// IWYU pragma: end_exports

// IWYU pragma: begin_keep
#include "pinocchio/algorithm/joint-configuration.hpp"
// IWYU pragma: end_keep
