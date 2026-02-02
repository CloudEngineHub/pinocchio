//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2015-2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_spatial_motion_hpp__
#define __pinocchio_spatial_motion_hpp__

// IWYU pragma: begin_keep
#include <Eigen/Core>

#include "pinocchio/common-traits.hpp"
#include "pinocchio/eigen-common.hpp"

#include "pinocchio/math/matrix.hpp"

#include "pinocchio/context.hxx"
#include "pinocchio/spatial/fwd.hxx"

#include "pinocchio/spatial/se3-common.hxx"
#include "pinocchio/spatial/skew.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/spatial/motion-common.hxx"
#include "pinocchio/spatial/motion-base.hxx"
#include "pinocchio/spatial/motion-dense.hxx"
#include "pinocchio/spatial/motion-tpl.hxx"
#include "pinocchio/spatial/motion-ref.hxx"
#include "pinocchio/spatial/motion-zero.hxx"
// IWYU pragma: end_exports

#endif // ifndef __pinocchio_spatial_motion_hpp__
