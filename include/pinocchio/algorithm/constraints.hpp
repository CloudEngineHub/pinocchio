//
// Copyright (c) INRIA 2026
//

#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <cstddef>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>
#include <utility>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <list>

#include <boost/blank.hpp>
#include <boost/core/ref.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/variant.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/contains.hpp>

#include <Eigen/Core>

#include "pinocchio/fwd.hpp"
#include "pinocchio/macros.hpp"
#include "pinocchio/unsupported.hpp"
#include "pinocchio/eigen-common.hpp"
#include "pinocchio/context.hpp"

#include "pinocchio/utils/check.hpp"
#include "pinocchio/utils/static-if.hpp"
#include "pinocchio/utils/reference.hpp"
#include "pinocchio/utils/std-vector.hpp"

#include "pinocchio/common.hpp"

#include "pinocchio/container/eigen-storage.hpp"
#include "pinocchio/container/matrix-stack.hpp"

#include "pinocchio/math.hpp"
#include "pinocchio/spatial.hpp"
#include "pinocchio/multibody.hpp"
#include "pinocchio/multibody/joint.hpp"
#include "pinocchio/serialization/serializable.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/check-model.hpp"
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

#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hxx"

#include "pinocchio/algorithm/constraints/blank-constraint.hxx"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hxx"

#include "pinocchio/algorithm/constraints/constraint-model-base.hxx"
#include "pinocchio/algorithm/constraints/constraint-data-base.hxx"

#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hxx"

#include "pinocchio/algorithm/constraints/constraint-data-generic.hxx"
#include "pinocchio/algorithm/constraints/constraint-model-generic.hxx"

#include "pinocchio/algorithm/constraints/kinematics-constraint-model-base.hxx"
#include "pinocchio/algorithm/constraints/binary-kinematics-constraint-model-base.hxx"

#include "pinocchio/algorithm/constraints/point-constraint-data-base.hxx"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hxx"
#include "pinocchio/algorithm/constraints/point-anchor-constraint.hxx"
#include "pinocchio/algorithm/constraints/point-contact-constraint.hxx"

#include "pinocchio/algorithm/constraints/frame-constraint-data-base.hxx"
#include "pinocchio/algorithm/constraints/frame-constraint-model-base.hxx"
#include "pinocchio/algorithm/constraints/frame-anchor-constraint.hxx"

#include "pinocchio/algorithm/constraints/jointwise-constraint-model-base.hxx"
#include "pinocchio/algorithm/constraints/joint-friction-constraint.hxx"
#include "pinocchio/algorithm/constraints/joint-limit-constraint.hxx"

#include "pinocchio/algorithm/constraints/constraint-collection-default.hxx"

#include "pinocchio/algorithm/constraints/constraint-ordering.hxx"
#include "pinocchio/algorithm/constraints/contact-info.hxx"

#include "pinocchio/algorithm/constraints/utils.hxx"
// IWYU pragma: end_exports
