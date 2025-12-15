//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_solvers_fwd_hpp__
#define __pinocchio_algorithm_solvers_fwd_hpp__

#include "pinocchio/algorithm/fwd.hpp"

namespace pinocchio
{
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
