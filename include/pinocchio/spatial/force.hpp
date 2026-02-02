//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_spatial_force_hpp__
#define __pinocchio_spatial_force_hpp__

// IWYU pragma: begin_keep
#include <Eigen/Core>

#include "pinocchio/common-traits.hpp"
#include "pinocchio/eigen-common.hpp"

#include "pinocchio/math/matrix.hpp"

#include "pinocchio/context.hxx" // IWYU pragma: keep
#include "pinocchio/spatial/fwd.hxx"
#include "pinocchio/spatial/se3-common.hxx"
#include "pinocchio/spatial/motion-common.hxx"
#include "pinocchio/spatial/force-common.hxx"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/spatial/force-common.hxx"
#include "pinocchio/spatial/force-base.hxx"
#include "pinocchio/spatial/force-dense.hxx"
#include "pinocchio/spatial/force-tpl.hxx"
#include "pinocchio/spatial/force-ref.hxx"
// IWYU pragma: end_exports

#endif // ifndef __pinocchio_spatial_force_hpp__
