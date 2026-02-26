//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_solvers_constraint_solver_base_hpp__
#define __pinocchio_algorithm_solvers_constraint_solver_base_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"

#ifdef PINOCCHIO_WITH_COLLISION
  #include <coal/timings.h>
#endif // PINOCCHIO_WITH_COLLISION

namespace pinocchio
{

  // fwd declaration for below
  template<typename Scalar>
  struct ConstraintSolverBaseTpl;

  template<typename _Scalar>
  struct ConstraintSolverResultBaseTpl;

  ///
  /// \brief Base struct for settings to pass to the solve method of a constraint solver.
  template<typename _Scalar>
  struct ConstraintSolverSettingsBaseTpl
  {
    typedef _Scalar Scalar;

    /// \brief Default constructor
    ConstraintSolverSettingsBaseTpl(
      std::size_t max_iterations,
      Scalar absolute_feasibility_tol,
      Scalar relative_feasibility_tol,
      Scalar absolute_complementarity_tol,
      Scalar relative_complementarity_tol,
      bool solve_ncp,
      bool measure_timings,
      bool stat_record)
    : max_iterations(max_iterations)
    , absolute_feasibility_tol(absolute_feasibility_tol)
    , relative_feasibility_tol(relative_feasibility_tol)
    , absolute_complementarity_tol(absolute_complementarity_tol)
    , relative_complementarity_tol(relative_complementarity_tol)
    , solve_ncp(solve_ncp)
    , measure_timings(measure_timings)
    , stat_record(stat_record)
    {
    }

    /// \brief Throws if settings are not valid.
    void checkValidity() const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        check_expression_if_real<Scalar>(absolute_feasibility_tol >= Scalar(0)),
        "absolute_feasibility_tol should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        check_expression_if_real<Scalar>(relative_feasibility_tol >= Scalar(0)),
        "relative_feasibility_tol should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        check_expression_if_real<Scalar>(absolute_complementarity_tol >= Scalar(0)),
        "absolute_complementarity_tol should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        check_expression_if_real<Scalar>(relative_complementarity_tol >= Scalar(0)),
        "relative_complementarity_tol should be >= 0.");
    }

    /// \brief Maximum number of iterations of the solver.
    std::size_t max_iterations;

    /// \brief Tolerance on the primal/dual feasibility.
    Scalar absolute_feasibility_tol;

    /// \brief Relative tolerance on the primal/dual feasibility.
    Scalar relative_feasibility_tol;

    /// \brief Absolute tolerance on the complementarity (duality gap).
    Scalar absolute_complementarity_tol;

    /// \brief Relative tolerance on the complementarity (duality gap).
    Scalar relative_complementarity_tol;

    /// \brief Whether or not to solve the NCP. If set to solve, the equivalent CCP
    /// is solved.
    bool solve_ncp;

    /// \brief Measure solve timings
    bool measure_timings;

    /// \brief Record per iteration stats.
    bool stat_record;
  }; // struct ConstraintSolverSettingsBaseTpl

  /// \brief Unsafe version of ConstraintSolverResultBaseTpl
  template<typename _Scalar>
  struct Unsafe<ConstraintSolverResultBaseTpl<_Scalar>>
  {
    typedef ConstraintSolverResultBaseTpl<_Scalar> SafeSelf;

    explicit Unsafe(SafeSelf & self)
    : self(self)
    {
    }

    /// \brief Make the solution valid.
    /// This is typically called by solvers when they are done with their computation.
    /// This allows them to mark the result as valid.
    void makeValid()
    {
      self.m_is_valid = true;
    }

  protected:
    SafeSelf & self;
  };

  ///
  /// \brief Base struct for result of a constraint solver.
  template<typename _Scalar>
  struct ConstraintSolverResultBaseTpl
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverResultBaseTpl Self;

    /// \brief Cast this class to its unsafe version.
    Unsafe<Self> unsafe()
    {
      return Unsafe<Self>(*this);
    }
    friend struct Unsafe<Self>;

    /// \brief Default constructor.
    ConstraintSolverResultBaseTpl()
    : iterations(0)
    , converged(false)
    , primal_feasibility(std::numeric_limits<Scalar>::quiet_NaN())
    , dual_feasibility(std::numeric_limits<Scalar>::quiet_NaN())
    , complementarity(std::numeric_limits<Scalar>::quiet_NaN())
    , m_is_valid(false)
    {
    }

    /// \brief Reset the solution.
    void reset()
    {
      iterations = 0;
      converged = false;
      primal_feasibility = std::numeric_limits<Scalar>::quiet_NaN();
      dual_feasibility = std::numeric_limits<Scalar>::quiet_NaN();
      complementarity = std::numeric_limits<Scalar>::quiet_NaN();

      m_is_valid = false;
    }

    /// \brief Returns true if the solution is in a valid state.
    /// If it is, it represents the result of a meaningful computation.
    /// Otherwise, it does not.
    bool isValid() const
    {
      return m_is_valid;
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

  protected:
    /// \brief Whether or not the solution is in a valid state.
    /// If it is, it represents the result of a meaningful computation.
    /// Otherwise, it does not.
    bool m_is_valid;

  }; // struct ConstraintSolverResultBaseTpl

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
      primal_feasibility.reserve(max_iterations);
      dual_feasibility.reserve(max_iterations);
      dual_feasibility_ncp.reserve(max_iterations);
      complementarity.reserve(max_iterations);
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

#ifdef PINOCCHIO_WITH_COLLISION
    typedef coal::CPUTimes CPUTimes;
    typedef coal::Timer Timer;
#endif // PINOCCHIO_WITH_COLLISION

    ConstraintSolverBaseTpl()
#ifdef PINOCCHIO_WITH_COLLISION
    : timer(false)
#endif // PINOCCHIO_WITH_COLLISION
    {
    }

#ifdef PINOCCHIO_WITH_COLLISION
    CPUTimes getCPUTimes() const
    {
      return timer.elapsed();
    }
#endif // PINOCCHIO_WITH_COLLISION

  protected:
#ifdef PINOCCHIO_WITH_COLLISION
    Timer timer;
#endif // PINOCCHIO_WITH_COLLISION

  }; // struct ConstraintSolverBaseTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_constraint_solver_base_hpp__
