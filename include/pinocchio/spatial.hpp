//
// Copyright (c) INRIA 2026
//

// IWYU pragma: begin_keep
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "pinocchio/traits.hpp"
#include "pinocchio/eigen-common.hpp"
#include "pinocchio/context.hxx"

#include "pinocchio/utils/cast.hpp"
#include "pinocchio/utils/static-if.hpp"
#include "pinocchio/utils/check.hpp"

#include "pinocchio/math.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/spatial/fwd.hxx"

#include "pinocchio/spatial/cartesian-axis.hxx"
#include "pinocchio/spatial/skew.hxx"
#include "pinocchio/spatial/symmetric3.hxx"

#include "pinocchio/spatial/se3-common.hxx"
#include "pinocchio/spatial/se3-base.hxx"
#include "pinocchio/spatial/se3-tpl.hxx"

#include "pinocchio/spatial/motion-common.hxx"
#include "pinocchio/spatial/motion-base.hxx"
#include "pinocchio/spatial/motion-dense.hxx"
#include "pinocchio/spatial/motion-tpl.hxx"
#include "pinocchio/spatial/motion-ref.hxx"
#include "pinocchio/spatial/motion-zero.hxx"

#include "pinocchio/spatial/spatial-axis.hxx"

#include "pinocchio/spatial/force-common.hxx"
#include "pinocchio/spatial/force-base.hxx"
#include "pinocchio/spatial/force-dense.hxx"
#include "pinocchio/spatial/force-tpl.hxx"
#include "pinocchio/spatial/force-ref.hxx"

#include "pinocchio/spatial/inertia.hxx"

#include "pinocchio/spatial/log-common.hxx"
#include "pinocchio/spatial/explog.hxx"
#include "pinocchio/spatial/explog-quaternion.hxx"
#include "pinocchio/spatial/log.hxx"
#include "pinocchio/spatial/se3-tpl-interpolate.hxx"

#include "pinocchio/spatial/classic-acceleration.hxx"

#include "pinocchio/spatial/act-on-set.hxx"
// IWYU pragma: end_exports
