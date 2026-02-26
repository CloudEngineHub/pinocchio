//
// Copyright (c) 2024 INRIA
//

#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <new>
#include <functional>
#include <utility>
#include <vector>
#include <type_traits>

#include <boost/fusion/container/vector.hpp>
#include <boost/container/detail/std_fwd.hpp>
#include <boost/variant.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "pinocchio/fwd.hpp"
#include "pinocchio/macros.hpp"
#include "pinocchio/eigen-common.hpp"

#include "pinocchio/utils/alloca.hpp"
#include "pinocchio/utils/promote-static-eval.hpp"
#include "pinocchio/utils/reference.hpp"
#include "pinocchio/utils/size-in-bytes.hpp"
#include "pinocchio/utils/template-template-parameter.hpp"

#include "pinocchio/container/eigen-storage.hpp"

#include "pinocchio/math.hpp"
#include "pinocchio/math/matrix-block.hpp"
#include "pinocchio/math/arithmetic-operators.hpp"
#include "pinocchio/multibody.hpp"
#include "pinocchio/multibody/joint.hpp"

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/preconditioner-base.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/algorithm/delassus-operator-base.hxx"
#include "pinocchio/algorithm/delassus-operator-sparse.hxx"
#include "pinocchio/algorithm/delassus-operator-rigid-body-visitors.hxx"
#include "pinocchio/algorithm/delassus-operator-rigid-body.hxx"
#include "pinocchio/algorithm/delassus-operator-dense.hxx"
#include "pinocchio/algorithm/delassus-operator-preconditioned.hxx"
#include "pinocchio/algorithm/delassus-operator-cholesky-expression.hxx"
// IWYU pragma: end_exports
