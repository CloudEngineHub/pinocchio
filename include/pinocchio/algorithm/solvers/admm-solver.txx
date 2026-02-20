//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_algorithm_solvers_admm_solver_txx__
#define __pinocchio_algorithm_solvers_admm_solver_txx__

#ifndef PINOCCHIO_SKIP_ALGORITHM_SOLVERS

  #include <functional>

namespace pinocchio
{

  // -------------------------------------------------------------------------
  // Struct instantiations
  // -------------------------------------------------------------------------

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    ADMMSolverSettingsTpl<context::Scalar>;

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    ADMMSolverStatsTpl<context::Scalar>;

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    ADMMSolverResultTpl<context::Scalar, context::Options>;

  namespace internal
  {
    extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
      ADMMSolverWorkspaceTpl<context::Scalar, context::Options>;
  } // namespace internal

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    ADMMConstraintSolverTpl<context::Scalar, context::Options>;

  // -------------------------------------------------------------------------
  // solve() with DelassusOperatorDense + default constraint collection
  // -------------------------------------------------------------------------

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI bool
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

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI bool
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

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI bool
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

#endif // ifndef PINOCCHIO_SKIP_ALGORITHM_SOLVERS

#endif // ifndef __pinocchio_algorithm_solvers_admm_solver_txx__
