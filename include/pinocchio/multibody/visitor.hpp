//
// Copyright (c) 2015-2019 CNRS INRIA
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#pragma once

// IWYU pragma: begin_keep
#include "pinocchio/traits.hpp"

#define BOOST_FUSION_INVOKE_MAX_ARITY 12
#include <boost/blank.hpp>
#include <boost/fusion/include/invoke.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/variant.hpp>

#include "pinocchio/context.hxx"

#include "pinocchio/math.hpp"

#include "pinocchio/multibody/fwd.hxx"
#include "pinocchio/multibody/joint/fwd.hxx"
#include "pinocchio/multibody/joint/joint-model-base.hxx"
#include "pinocchio/multibody/joint/joint-data-base.hxx"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/multibody/visitor/fusion.hxx"
#include "pinocchio/multibody/visitor/joint-unary-visitor.hxx"
#include "pinocchio/multibody/visitor/joint-binary-visitor.hxx"
// IWYU pragma: end_exports
