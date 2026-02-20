//
// Copyright (c) 2026 INRIA
//

#include "pinocchio/algorithm/solvers/clarabel-solver.hpp"

#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT

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
    ClarabelSolverSettingsTpl<context::Scalar>;

  template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    ClarabelSolverResultTpl<context::Scalar, context::Options>;

  template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
    ClarabelConstraintSolverTpl<context::Scalar, context::Options>;

  namespace internal
  {
    template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI
      ClarabelSolverWorkspaceTpl<context::Scalar, context::Options>;
  } // namespace internal

  // -------------------------------------------------------------------------
  // solve() with DelassusOperatorDense + default constraint collection
  // -------------------------------------------------------------------------

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  ClarabelConstraintSolverTpl<context::Scalar, context::Options>::solve<
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
    const ClarabelSolverSettingsTpl<context::Scalar> &,
    ClarabelSolverResultTpl<context::Scalar, context::Options> &);

  // -------------------------------------------------------------------------
  // solve() with DelassusCholeskyExpression + default constraint collection
  // -------------------------------------------------------------------------

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  ClarabelConstraintSolverTpl<context::Scalar, context::Options>::solve<
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
    const ClarabelSolverSettingsTpl<context::Scalar> &,
    ClarabelSolverResultTpl<context::Scalar, context::Options> &);

  // -------------------------------------------------------------------------
  // solve() with DelassusOperatorRigidBodySystems + default constraint collection
  // -------------------------------------------------------------------------

  template PINOCCHIO_EXPLICIT_INSTANTIATION_DEFINITION_DLLAPI bool
  ClarabelConstraintSolverTpl<context::Scalar, context::Options>::solve<
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
    const ClarabelSolverSettingsTpl<context::Scalar> &,
    ClarabelSolverResultTpl<context::Scalar, context::Options> &);

} // namespace pinocchio

#endif // PINOCCHIO_WITH_CLARABEL_SUPPORT
