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

#include <boost/type_traits.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/number.hpp>

#include "pinocchio/context.hxx"
#include "pinocchio/multibody/fwd.hpp"

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

#include "pinocchio/math/matrix-inverse-code-generated.hxx"
#include "pinocchio/math/details/matrix-inverse-1x1.hxx"
#include "pinocchio/math/details/matrix-inverse-2x2.hxx"
#include "pinocchio/math/details/matrix-inverse-3x3.hxx"
#include "pinocchio/math/details/matrix-inverse-4x4.hxx"
#include "pinocchio/math/details/matrix-inverse-5x5.hxx"
#include "pinocchio/math/details/matrix-inverse-6x6.hxx"
#include "pinocchio/math/details/matrix-inverse-7x7.hxx"
#include "pinocchio/math/details/matrix-inverse-8x8.hxx"
#include "pinocchio/math/details/matrix-inverse-9x9.hxx"
#include "pinocchio/math/details/matrix-inverse-10x10.hxx"
#include "pinocchio/math/details/matrix-inverse-11x11.hxx"
#include "pinocchio/math/details/matrix-inverse-12x12.hxx"
#include "pinocchio/math/matrix-inverse.hxx"
#include "pinocchio/math/matrix-product.hxx"
// IWYU pragma: end_exports
