//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <limits>
#include <ostream>
#include <string>
#include <vector>
#include <type_traits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#define BOOST_FUSION_INVOKE_MAX_ARITY 12
#include <boost/type_traits.hpp>
#include <boost/blank.hpp>
#include <boost/fusion/include/invoke.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/get.hpp>

#include "pinocchio/context.hxx"
#include "pinocchio/traits.hpp"
#include "pinocchio/macros.hpp"

#include "pinocchio/utils/static-if.hpp"
#include "pinocchio/utils/axis-label.hpp"

#include "pinocchio/math.hpp"
#include "pinocchio/spatial.hpp"

#include "pinocchio/multibody/fwd.hxx"
#include "pinocchio/multibody/joint-motion-subspace-base.hxx"
#include "pinocchio/multibody/joint-motion-subspace-generic.hxx"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/multibody/joint/fwd.hxx"
#include "pinocchio/multibody/joint/joint-model-base.hxx"
#include "pinocchio/multibody/joint/joint-data-base.hxx"
#include "pinocchio/multibody/joint/joint-common-operations.hxx"

#include "pinocchio/multibody/visitor/fusion.hxx"
#include "pinocchio/multibody/visitor/joint-unary-visitor.hxx"
#include "pinocchio/multibody/visitor/joint-binary-visitor.hxx"

#include "pinocchio/multibody/joint/joint-revolute.hxx"
#include "pinocchio/multibody/joint/joint-revolute-unaligned.hxx"
#include "pinocchio/multibody/joint/joint-revolute-unbounded.hxx"
#include "pinocchio/multibody/joint/joint-revolute-unbounded-unaligned.hxx"

#include "pinocchio/multibody/joint/joint-translation.hxx"

#include "pinocchio/multibody/joint/joint-prismatic.hxx"
#include "pinocchio/multibody/joint/joint-prismatic-unaligned.hxx"

#include "pinocchio/multibody/joint/joint-planar.hxx"

#include "pinocchio/multibody/joint/joint-ellipsoid.hxx"

#include "pinocchio/multibody/joint/joint-helical.hxx"
#include "pinocchio/multibody/joint/joint-helical-unaligned.hxx"

#include "pinocchio/multibody/joint/joint-spherical.hxx"
#include "pinocchio/multibody/joint/joint-spherical-ZYX.hxx"

#include "pinocchio/multibody/joint/joint-universal.hxx"

#include "pinocchio/multibody/joint/joint-free-flyer.hxx"
// IWYU pragma: end_exports
