//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_solvers_constraint_solver_base_hpp__
#define __pinocchio_algorithm_solvers_constraint_solver_base_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"

#ifdef PINOCCHIO_WITH_HPP_FCL
  #include <hpp/fcl/timings.h>
#endif // PINOCCHIO_WITH_HPP_FCL

namespace pinocchio
{

  ///
  /// \brief Base struct for settings to pass to the solve method of a constraint solver.
  template<typename _Scalar>
  struct ConstraintSolverSettingsBaseTpl
  {
    typedef _Scalar Scalar;

    /// \brief Default constructor
    ConstraintSolverSettingsBaseTpl(
      std::size_t max_iterations,
      Scalar absolute_tol_feasibility,
      Scalar relative_tol_feasibility,
      Scalar absolute_tol_complementarity,
      Scalar relative_tol_complementarity,
      bool solve_ncp,
      bool measure_timings,
      bool stat_record)
    : max_iterations(max_iterations)
    , absolute_tol_feasibility(absolute_tol_feasibility)
    , relative_tol_feasibility(relative_tol_feasibility)
    , absolute_tol_complementarity(absolute_tol_complementarity)
    , relative_tol_complementarity(relative_tol_complementarity)
    , solve_ncp(solve_ncp)
    , measure_timings(measure_timings)
    , stat_record(stat_record)
    {
    }

    /// \brief Throws if settings are not valid.
    void checkValidity() const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        absolute_tol_feasibility >= Scalar(0), "absolute_tol_feasibility should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        relative_tol_feasibility >= Scalar(0), "relative_tol_feasibility should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        absolute_tol_complementarity >= Scalar(0), "absolute_tol_complementarity should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        relative_tol_complementarity >= Scalar(0), "relative_tol_complementarity should be >= 0.");
    }

    /// \brief Maximum number of iterations of the solver.
    std::size_t max_iterations;

    /// \brief Tolerance on the primal/dual feasibility.
    Scalar absolute_tol_feasibility;

    /// \brief Relative tolerance on the primal/dual feasibility.
    Scalar relative_tol_feasibility;

    /// \brief Absolute tolerance on the complementarity (duality gap).
    Scalar absolute_tol_complementarity;

    /// \brief Relative tolerance on the complementarity (duality gap).
    Scalar relative_tol_complementarity;

    /// \brief Whether or not to solve the NCP. If set to solve, the equivalent CCP
    /// is solved.
    bool solve_ncp;

    /// \brief Measure solve timings
    bool measure_timings;

    /// \brief Record per iteration stats.
    bool stat_record;
  }; // struct ConstraintSolverSettingsBaseTpl

  ///
  /// \brief Base struct for settings to pass to the solve method of a constraint solver.
  template<typename _Scalar>
  struct ConstraintSolverSolutionBaseTpl
  {
    typedef _Scalar Scalar;
    static constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();

    /// \brief Default constructor.
    ConstraintSolverSolutionBaseTpl()
    : iterations(0)
    , converged(false)
    , primal_feasibility(nan)
    , dual_feasibility(nan)
    , complementarity(nan)
    {
    }

    /// \brief Reset the solution.
    void reset()
    {
      iterations = 0;
      converged = false;
      primal_feasibility = nan;
      dual_feasibility = nan;
      complementarity = nan;
    }

    /// \brief Number of iterations of the solver
    std::size_t iterations;

    /// \brief Whether or not the solver has converged
    bool converged;

    /// \brief Value of the primal feasibility
    Scalar primal_feasibility;

    /// \brief Value of the dual feasibility
    Scalar dual_feasibility;

    /// \brief Value of the complementarity
    Scalar complementarity;

  }; // struct ConstraintSolverSolutionBaseTpl

  ///
  /// \brief Base struct to track a constraint solver progress per iteration.
  template<typename _Scalar>
  struct ConstraintSolverStatsBaseTpl
  {
    typedef _Scalar Scalar;

    /// \brief Default constructor.
    ConstraintSolverStatsBaseTpl()
    : iterations(0)
    {
    }

    /// \brief Constructor given a maximum iteration of the solver.
    explicit ConstraintSolverStatsBaseTpl(std::size_t max_iterations)
    : iterations(0)
    {
      reserve(max_iterations);
    }

    /// \brief Reserve enough storage for max_it iterations.
    void reserve(std::size_t max_iterations)
    {
      primal_feasibility.reserve(size_t(max_iterations));
      dual_feasibility.reserve(size_t(max_iterations));
      dual_feasibility_ncp.reserve(size_t(max_iterations));
      complementarity.reserve(size_t(max_iterations));
    }

    /// \brief Reset stats.
    void reset()
    {
      iterations = 0;
      primal_feasibility.clear();
      dual_feasibility.clear();
      dual_feasibility_ncp.clear();
      complementarity.clear();
    }

    /// \brief Returns the size of the stats (number of iterations tracked).
    size_t size() const
    {
      return primal_feasibility.size();
    }

    ///  \brief Total number of iterations.
    std::size_t iterations;

    /// \brief History of primal feasibility values.
    std::vector<Scalar> primal_feasibility;

    /// \brief History of dual feasibility values.
    std::vector<Scalar> dual_feasibility;

    /// \brief History of NCP dual feasibility values.
    std::vector<Scalar> dual_feasibility_ncp;

    /// \brief History of complementarity values.
    std::vector<Scalar> complementarity;
  };

  template<typename _Scalar>
  struct ConstraintSolverBaseTpl
  {
    typedef _Scalar Scalar;

#ifdef PINOCCHIO_WITH_HPP_FCL
    typedef hpp::fcl::CPUTimes CPUTimes;
    typedef hpp::fcl::Timer Timer;
#endif // PINOCCHIO_WITH_HPP_FCL

    ConstraintSolverBaseTpl()
    :
#ifdef PINOCCHIO_WITH_HPP_FCL
      timer(false)
#endif // PINOCCHIO_WITH_HPP_FCL
    {
    }

#ifdef PINOCCHIO_WITH_HPP_FCL
    CPUTimes getCPUTimes() const
    {
      return timer.elapsed();
    }
#endif // PINOCCHIO_WITH_HPP_FCL

  protected:
#ifdef PINOCCHIO_WITH_HPP_FCL
    Timer timer;
#endif // PINOCCHIO_WITH_HPP_FCL

  }; // struct ConstraintSolverBaseTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_constraint_solver_base_hpp__
