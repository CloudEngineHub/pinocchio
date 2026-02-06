//
// Copyright (c) 2026 INRIA
//

#pragma once

#define PINOCCHIO_SKIP_ALGORITHM_MODEL
#define PINOCCHIO_SKIP_ALGORITHM_GEOMETRY
#define PINOCCHIO_SKIP_MULTIBODY_SAMPLE_MODELS
#define PINOCCHIO_TEMPLATE_INSTANTIATION_HEADER "pinocchio/codegen/cppadcg.hpp"

#include <cppad/cg.hpp>
#define PINOCCHIO_SCALAR_TYPE ::CppAD::AD<CppAD::cg::CG<double>>
#include "pinocchio/context/generic.hxx"
#undef PINOCCHIO_SCALAR_TYPE
