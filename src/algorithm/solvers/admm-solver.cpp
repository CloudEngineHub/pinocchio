//
// Copyright (c) 2026 INRIA
//

#include "pinocchio/algorithm/solvers/admm-solver.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include "pinocchio/algorithm/constraint-cholesky.hpp"
#include "pinocchio/algorithm/delassus-operator-cholesky-expression.hpp"
#include "pinocchio/algorithm/delassus-operator-rigid-body.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"

namespace pinocchio
{

  // -------------------------------------------------------------------------
  // Struct instantiations
  // -------------------------------------------------------------------------

  template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    ADMMSolverSettingsTpl<context::Scalar>;

  template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    ADMMSolverStatsTpl<context::Scalar>;

  template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    ADMMSolverResultTpl<context::Scalar, context::Options>;

  namespace internal
  {
    template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
      ADMMSolverWorkspaceTpl<context::Scalar, context::Options>;
  } // namespace internal

  template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    ADMMConstraintSolverTpl<context::Scalar, context::Options>;

  // -------------------------------------------------------------------------
  // solve() with DelassusOperatorDense + default constraint collection
  // -------------------------------------------------------------------------

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  ADMMConstraintSolverTpl<context::Scalar, context::Options>::solve<
    DelassusOperatorDenseTpl<context::Scalar, context::Options>,
    context::VectorXs,
    ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>,
    ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>>(
    DelassusOperatorBase<DelassusOperatorDenseTpl<context::Scalar, context::Options>> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const std::vector<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const std::vector<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const ADMMSolverSettingsTpl<context::Scalar> &,
    ADMMSolverResultTpl<context::Scalar, context::Options> &);

  // -------------------------------------------------------------------------
  // solve() with DelassusCholeskyExpression + default constraint collection
  // -------------------------------------------------------------------------

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  ADMMConstraintSolverTpl<context::Scalar, context::Options>::solve<
    DelassusCholeskyExpressionTpl<
      ContactCholeskyDecompositionTpl<context::Scalar, context::Options>>,
    context::VectorXs,
    ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>,
    ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>>(
    DelassusOperatorBase<DelassusCholeskyExpressionTpl<
      ContactCholeskyDecompositionTpl<context::Scalar, context::Options>>> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const std::vector<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const std::vector<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const ADMMSolverSettingsTpl<context::Scalar> &,
    ADMMSolverResultTpl<context::Scalar, context::Options> &);

  // -------------------------------------------------------------------------
  // solve() with DelassusOperatorRigidBodySystems + default constraint collection
  // -------------------------------------------------------------------------

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  ADMMConstraintSolverTpl<context::Scalar, context::Options>::solve<
    DelassusOperatorRigidBodySystemsTpl<
      context::Scalar,
      context::Options,
      JointCollectionDefaultTpl,
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::reference_wrapper>,
    context::VectorXs,
    ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>,
    ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>>(
    DelassusOperatorBase<DelassusOperatorRigidBodySystemsTpl<
      context::Scalar,
      context::Options,
      JointCollectionDefaultTpl,
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::reference_wrapper>> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const std::vector<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const std::vector<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const ADMMSolverSettingsTpl<context::Scalar> &,
    ADMMSolverResultTpl<context::Scalar, context::Options> &);

} // namespace pinocchio
