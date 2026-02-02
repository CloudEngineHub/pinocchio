//
// Copyright (c) 2015-2020 CNRS
//

#ifndef __pinocchio_spatial_se3_hpp__
#define __pinocchio_spatial_se3_hpp__

// IWYU pragma: begin_keep
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "pinocchio/eigen-common.hpp"
#include "pinocchio/common-traits.hpp"

// TODO replace by fwd.hpp
#include "pinocchio/context.hxx" // IWYU pragma: keep
#include "pinocchio/spatial/fwd.hxx"

#include "pinocchio/utils/cast.hpp"

#include "pinocchio/math/matrix.hpp"
#include "pinocchio/math/quaternion.hpp"
#include "pinocchio/math/rotation.hpp"

#include "pinocchio/spatial/cartesian-axis.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/spatial/se3-common.hxx"
#include "pinocchio/spatial/se3-base.hxx"
#include "pinocchio/spatial/se3-tpl.hxx"
// IWYU pragma: end_exports

#endif // ifndef __pinocchio_spatial_se3_hpp__
