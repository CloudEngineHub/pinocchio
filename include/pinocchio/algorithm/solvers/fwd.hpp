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
  template<typename Scalar>
  struct PGSConstraintSolverTpl;
  typedef PGSConstraintSolverTpl<context::Scalar> PGSConstraintSolver;

  template<typename Scalar>
  struct PGSSolverSettingsTpl;
  typedef PGSSolverSettingsTpl<context::Scalar> PGSSolverSettings;

  template<typename Scalar>
  struct PGSSolverSolutionTpl;
  typedef PGSSolverSolutionTpl<context::Scalar> PGSSolverSolution;

  template<typename Scalar>
  struct PGSSolverStatsTpl;
  typedef PGSSolverStatsTpl<context::Scalar> PGSSolverStats;

  // -------------------------
  // ADMM constraint solver
  template<typename Scalar>
  struct ADMMConstraintSolverTpl;
  typedef ADMMConstraintSolverTpl<context::Scalar> ADMMConstraintSolver;

  template<typename Scalar>
  struct ADMMSolverSettingsTpl;
  typedef ADMMSolverSettingsTpl<context::Scalar> ADMMSolverSettings;

  template<typename Scalar>
  struct ADMMSolverSolutionTpl;
  typedef ADMMSolverSolutionTpl<context::Scalar> ADMMSolverSolution;

  template<typename Scalar>
  struct ADMMSolverStatsTpl;
  typedef ADMMSolverStatsTpl<context::Scalar> ADMMSolverStats;

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_fwd_hpp__
