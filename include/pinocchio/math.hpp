//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <unsupported/Eigen/CXX11/Tensor>

#include <cassert>
#include <cmath>
#include <limits>
#include <cstddef>
#include <type_traits>

#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/detail/number_base.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/fwd.hpp>
#include <boost/random.hpp>

#include "pinocchio/context.hxx"
#include "pinocchio/multibody/fwd.hxx"

#include "pinocchio/macros.hpp"
#include "pinocchio/eigen-common.hpp"
#include "pinocchio/alloca.hpp"

#include "pinocchio/utils/static-if.hpp"
#include "pinocchio/utils/check.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/math/comparison-operators.hxx"
#include "pinocchio/math/fwd.hxx"
#include "pinocchio/math/alias.hxx"
#include "pinocchio/math/taylor-expansion.hxx"
#include "pinocchio/math/matrix.hxx"
#include "pinocchio/math/sincos.hxx"
#include "pinocchio/math/rpy.hxx"
#include "pinocchio/math/rotation.hxx"
#include "pinocchio/math/quaternion.hxx"
#include "pinocchio/math/triangular-matrix.hxx"
#include "pinocchio/math/matrix-block.hxx"
#include "pinocchio/math/tensor.hxx"
#include "pinocchio/math/eigenvalues.hxx"
#include "pinocchio/math/multiprecision.hxx"
#include "pinocchio/math/gram-schmidt-orthonormalisation.hxx"
#include "pinocchio/math/sign.hxx"
#include "pinocchio/math/eigenvalues-tridiagonal-matrix.hxx"
#include "pinocchio/math/tridiagonal-matrix.hxx"
#include "pinocchio/math/lanczos-decomposition.hxx"

// IWYU pragma: end_exports
