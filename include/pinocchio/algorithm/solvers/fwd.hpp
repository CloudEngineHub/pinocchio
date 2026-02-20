//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_solvers_fwd_hpp__
#define __pinocchio_algorithm_solvers_fwd_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{
  // -------------------------
  // PGS constraint solver
  template<typename Scalar, int Options>
  struct PGSConstraintSolverTpl;
  typedef PGSConstraintSolverTpl<context::Scalar, context::Options> PGSConstraintSolver;

  template<typename Scalar>
  struct PGSSolverSettingsTpl;
  typedef PGSSolverSettingsTpl<context::Scalar> PGSSolverSettings;

  template<typename Scalar, int Options>
  struct PGSSolverResultTpl;
  typedef PGSSolverResultTpl<context::Scalar, context::Options> PGSSolverResult;

  template<typename Scalar>
  struct PGSSolverStatsTpl;
  typedef PGSSolverStatsTpl<context::Scalar> PGSSolverStats;

  // -------------------------
  // ADMM constraint solver
  template<typename Scalar, int Options>
  struct ADMMConstraintSolverTpl;
  typedef ADMMConstraintSolverTpl<context::Scalar, context::Options> ADMMConstraintSolver;

  template<typename Scalar>
  struct ADMMSolverSettingsTpl;
  typedef ADMMSolverSettingsTpl<context::Scalar> ADMMSolverSettings;

  template<typename Scalar, int Options>
  struct ADMMSolverResultTpl;
  typedef ADMMSolverResultTpl<context::Scalar, context::Options> ADMMSolverResult;

  template<typename Scalar>
  struct ADMMSolverStatsTpl;
  typedef ADMMSolverStatsTpl<context::Scalar> ADMMSolverStats;

#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT
  // -------------------------
  // Optional clarabel IPM constraint solver
  template<typename Scalar, int Options>
  struct ClarabelConstraintSolverTpl;
  typedef ClarabelConstraintSolverTpl<context::Scalar, context::Options> ClarabelConstraintSolver;

  template<typename Scalar>
  struct ClarabelSolverSettingsTpl;
  typedef ClarabelSolverSettingsTpl<context::Scalar> ClarabelSolverSettings;

  template<typename Scalar, int Options>
  struct ClarabelSolverResultTpl;
  typedef ClarabelSolverResultTpl<context::Scalar, context::Options> ClarabelSolverResult;
#endif

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_fwd_hpp__
