//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_algorithm_solvers_pgs_solver_txx__
#define __pinocchio_algorithm_solvers_pgs_solver_txx__

namespace pinocchio
{

  // -------------------------------------------------------------------------
  // Struct instantiations
  // -------------------------------------------------------------------------

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    PGSSolverSettingsTpl<context::Scalar>;

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    PGSSolverStatsTpl<context::Scalar>;

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    PGSSolverResultTpl<context::Scalar, context::Options>;

  namespace internal
  {
    extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
      PGSSolverWorkspaceTpl<context::Scalar, context::Options>;
  } // namespace internal

  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    PGSConstraintSolverTpl<context::Scalar, context::Options>;

  // -------------------------------------------------------------------------
  // solve() with MatrixXs + default constraint collection
  // -------------------------------------------------------------------------

  extern template PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI bool
  PGSConstraintSolverTpl<context::Scalar, context::Options>::solve<
    context::MatrixXs,
    context::VectorXs,
    ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>,
    ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
    std::allocator<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>>(
    const Eigen::MatrixBase<context::MatrixXs> &,
    const Eigen::MatrixBase<context::VectorXs> &,
    const std::vector<
      ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const std::vector<
      ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>,
      std::allocator<
        ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>>> &,
    const PGSSolverSettingsTpl<context::Scalar> &,
    PGSSolverResultTpl<context::Scalar, context::Options> &);

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_pgs_solver_txx__
