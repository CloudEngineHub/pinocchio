//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_solvers_pgs_solver_hpp__
#define __pinocchio_algorithm_solvers_pgs_solver_hpp__

#include "pinocchio/algorithm/solvers/fwd.hpp"
#include "pinocchio/algorithm/solvers/constraint-solver-base.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"

#include "pinocchio/container/eigen-storage.hpp"

#include <optional>
#include <limits>

namespace pinocchio
{
  // fwd declarations for PGS-internal structs
  // user-api structs are fwd delclared in solvers/fwd.hpp.
  // see below for definitions
  namespace internal
  {
    template<typename Scalar>
    struct PGSSolverWorkspaceTpl;
  }

  /// \brief Projected Gauss Siedel solver
  template<typename _Scalar>
  struct PGSConstraintSolverTpl : ConstraintSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;

    typedef internal::PGSSolverWorkspaceTpl<Scalar> PGSSolverWorkspace;
    typedef PGSSolverSettingsTpl<Scalar> PGSSolverSettings;
    typedef PGSSolverResultTpl<Scalar> PGSSolverResult;
    typedef PGSSolverStatsTpl<Scalar> PGSSolverStats;

    explicit PGSConstraintSolverTpl(std::size_t problem_size = 0)
    : Base()
    , stats()
    , m_workspace(problem_size)
    , m_is_valid(false)
    {
    }

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_models,
    /// constraint_datas).
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velocity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in] constraint_datas Vector of constraint datas.
    /// \param[in] settings Settings for the PGS solver.
    /// \param[in/out] result Solution to the constraint problem. Also contains the warmstart to
    /// solve the problem.
    ///
    /// \returns True if the problem has converged.
    template<
      typename MatrixType,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    bool solve(
      const Eigen::MatrixBase<MatrixType> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const PGSSolverSettings & settings,
      PGSSolverResult & result);

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_models,
    /// constraint_datas).
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velocity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in] constraint_datas Vector of constraint datas.
    /// \param[in] settings Settings for the PGS solver.
    /// \param[in/out] result Solution to the constraint problem. Also contains the warmstart to
    /// solve the problem.
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    bool solve(
      const DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const PGSSolverSettings & settings,
      PGSSolverResult & result)
    {
      return solve(
        delassus.derived().matrix(), g, constraint_models, constraint_datas, settings, result);
    }

    /// \brief Reset the constraint solver as if it has never run.
    void reset()
    {
      stats.reset();
      m_workspace.reset();
      m_is_valid = false;
    }

    /// \brief Returns true if solver is in a valid state (it has solved a constraint problem).
    /// If so, its stats are valid.
    bool isValid() const
    {
      return m_is_valid;
    }

#ifdef PINOCCHIO_WITH_HPP_FCL
    /// \brief Timer for the `solve` method
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL

    /// \brief Per-iteration stats of the PGS solver.
    PGSSolverStats stats;

  protected:
    /// \brief Workspace of the PGS solver.
    /// This is an internal of the solver and is not meant to be accessed by
    /// users.
    PGSSolverWorkspace m_workspace;

    /// \brief Flag to check whether or not the solver is in a reset state.
    /// If not, its stats are valid.
    bool m_is_valid;

  }; // struct PGSConstraintSolverTpl

  ///
  /// \brief Settings for the PGS constraint solver loop.
  template<typename _Scalar>
  struct PGSSolverSettingsTpl : ConstraintSolverSettingsBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverSettingsBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;

    /// \brief Default constructor
    PGSSolverSettingsTpl(
      std::size_t max_iterations = 1000,
      Scalar absolute_feasibility_tol = Scalar(1e-6),
      Scalar relative_feasibility_tol = Scalar(1e-6),
      Scalar absolute_complementarity_tol = Scalar(1e-6),
      Scalar relative_complementarity_tol = Scalar(1e-6),
      bool solve_ncp = true,
      bool measure_timings = false,
      bool stat_record = false,
      Scalar over_relaxation = Scalar(1))
    : Base(
        max_iterations,
        absolute_feasibility_tol,
        relative_feasibility_tol,
        absolute_complementarity_tol,
        relative_complementarity_tol,
        solve_ncp,
        measure_timings,
        stat_record)
    , over_relaxation(over_relaxation)
    {
    }

    /// \brief Throws if settings are not valid.
    void checkValidity() const
    {
      Base::checkValidity();
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        over_relaxation < Scalar(2) && over_relaxation > Scalar(0),
        "over_relaxation should lie in ]0,2[.");
    }

    // ----------------------
    // General settings

    /// \brief Maximum number of iterations of the solver.
    using Base::max_iterations;

    /// \brief Tolerance on the primal/dual feasibilityibility.
    using Base::absolute_feasibility_tol;

    /// \brief Relative tolerance on the primal/dual feasibilityibility.
    using Base::relative_feasibility_tol;

    /// \brief Absolute tolerance of complementarity (duality complementarity).
    using Base::absolute_complementarity_tol;

    /// \brief Relative tolerance of complementarity (duality complementarity).
    using Base::relative_complementarity_tol;

    /// \brief Whether or not to solve the NCP. If set to solve, the equivalent CCP
    /// is solved.
    using Base::solve_ncp;

    /// \brief Measure solve timings
    using Base::measure_timings;

    /// \brief Record per iteration stats.
    using Base::stat_record;

    // ----------------------
    // PGS specific settings

    /// \brief Over-relaxation of PGS step. Default value is 1.
    Scalar over_relaxation;
  }; // struct PGSSolverSettingsTpl

  ///
  /// \brief Struct describing the solution of the PGS constraint solver
  /// after calling the `solve` method.
  /// Also contains the warmstart of the solution to the constraint problem.
  template<typename _Scalar>
  struct PGSSolverResultTpl : ConstraintSolverResultBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverResultBaseTpl<Scalar> Base;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;
    typedef EigenStorageTpl<VectorXs> VectorXsStorage;

    // make PGS solver a friend so that it can use `makeValid`
    friend struct PGSConstraintSolverTpl<Scalar>;

    using Base::isValid;

    /// \brief Default constructor.
    PGSSolverResultTpl()
    : Base()
    , problem_size(0)
    , primal_guess(std::nullopt)
    , x(x_storage.map())
    , y(y_storage.map())
    {
    }

    /// \brief Reset the results.
    /// \note This method does not touch the warmstart fields.
    void reset(std::size_t problem_size_ = 0)
    {
      Base::reset();
      problem_size = problem_size_;

      primal_guess.reset();

      resize(problem_size);

      // set solution to nan - solver has not run
      x.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      y.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
    }

    /// \brief Resize the primal/dual vectors of the solution.
    void resize(std::size_t problem_size_)
    {
      problem_size = problem_size_;

      Eigen::Index np = static_cast<Eigen::Index>(problem_size);
      x_storage.resize(np);
      y_storage.resize(np);
    }

    /// \brief Retrieve primal solution.
    template<typename VectorLike>
    void retrievePrimalSolution(const Eigen::MatrixBase<VectorLike> & primal_solution_) const
    {
      auto & primal_solution = primal_solution_.const_cast_derived();
      primal_solution = x;
    }

    /// \brief Retrieve dual solution.
    template<typename VectorLike>
    void retrieveDualSolution(const Eigen::MatrixBase<VectorLike> & dual_solution_) const
    {
      auto & dual_solution = dual_solution_.const_cast_derived();
      dual_solution = y;
    }

    /// \brief Retrieve constraint impulses.
    template<typename VectorLike>
    void
    retrieveConstraintImpulses(const Eigen::MatrixBase<VectorLike> & constraint_impulses_) const
    {
      auto & constraint_impulses = constraint_impulses_.const_cast_derived();
      constraint_impulses = x;
    }

    /// \brief Retrieve constraint velocities.
    /// At the optimum we have y = Gx + g.
    /// WARNING: the PGS solver does not take into account desaxce terms for now.
    /// It only solves the CCP (not the NCP).
    template<typename VectorLike>
    void
    retrieveConstraintVelocities(const Eigen::MatrixBase<VectorLike> & constraint_velocities_) const
    {
      auto & constraint_velocities = constraint_velocities_.const_cast_derived();
      constraint_velocities = y;
    }

    /// \brief Number of iterations of the solver.
    using Base::iterations;

    /// \brief Whether or not the solver has converged.
    using Base::converged;

    /// \brief Value of the primal feasibility.
    using Base::primal_feasibility;

    /// \brief Value of the dual feasibility.
    using Base::dual_feasibility;

    /// \brief Value of the complementarity.
    using Base::complementarity;

    /// \brief Size of primal/dual variables.
    std::size_t problem_size;

    // ----------------------
    // Solution warmstart

    /// \brief Optional guess for the primal variable (impulses).
    std::optional<RefConstVectorXs> primal_guess;

    // ----------------------
    // Solution - output of the solver

    /// \brief Primal solution.
    /// \note Order of storage/map declaration is important!
    /// First declare the storage, then the map, otherwise map will point to nothing.
    VectorXsStorage x_storage;
    typename VectorXsStorage::RefMapType x;

    /// \brief Dual solution.
    VectorXsStorage y_storage;
    typename VectorXsStorage::RefMapType y;
  }; // struct PGSSolverResultTpl

  ///
  /// \brief Struct to track per iteration progress of PGS constraint solver.
  template<typename _Scalar>
  struct PGSSolverStatsTpl : ConstraintSolverStatsBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverStatsBaseTpl<Scalar> Base;

    /// \brief Default constructor.
    PGSSolverStatsTpl()
    : Base()
    {
    }

    /// \brief Constructor given a maximum iteration of the solver.
    explicit PGSSolverStatsTpl(std::size_t max_iterations)
    : Base(max_iterations)
    {
      Base::reserve(max_iterations);
    }

    /// \brief How many iterations the solver ran.
    using Base::iterations;

    /// \brief Vector storing per iteration primal feasibility.
    using Base::primal_feasibility;

    /// \brief Vector storing per iteration dual feasibility.
    using Base::dual_feasibility;

    /// \brief Vector storing per iteration NCP feasibility.
    using Base::dual_feasibility_ncp;

    /// \brief Vector storing per iteration complementarity.
    using Base::complementarity;
  }; // struct PGSSolverStatsTpl

  namespace internal
  {
    ///
    /// \brief Workspace for the PGS constraint solver.
    template<typename _Scalar>
    struct PGSSolverWorkspaceTpl
    {
      typedef _Scalar Scalar;
      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
      typedef EigenStorageTpl<VectorXs> VectorXsStorage;

      /// \brief Constructor given problem_size.
      PGSSolverWorkspaceTpl(std::size_t problem_size = 0)
      : problem_size(problem_size)
      , x(x_storage.map())
      , x_previous(x_previous_storage.map())
      , y(y_storage.map())
      , rhs(rhs_storage.map())
      , tmp(tmp_storage.map())
      {
        reset(problem_size);
      }

      /// \brief Reset the workspace.
      void reset(std::size_t problem_size_ = 0)
      {
        problem_size = problem_size_;

        resize(problem_size);

#ifndef NDEBUG
        x.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        x_previous.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        y.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        tmp.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        rhs.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
#endif
      }

      /// \brief Resize workspace vectors to problem size.
      void resize(std::size_t problem_size_)
      {
        problem_size = problem_size_;

        Eigen::Index np = static_cast<Eigen::Index>(problem_size);
        x_storage.resize(np);
        x_previous_storage.resize(np);
        y_storage.resize(np);
        tmp_storage.resize(np);
        rhs_storage.resize(np);
      }

      /// \brief Size of problem.
      std::size_t problem_size;

      /// \brief Primal variable (impulses) at current iteration.
      VectorXsStorage x_storage;
      typename VectorXsStorage::RefMapType x;

      /// \brief Primal variable (impulses) at previous iteration.
      VectorXsStorage x_previous_storage;
      typename VectorXsStorage::RefMapType x_previous;

      /// \brief Dual variable (constraint velocities) at current iteration.
      VectorXsStorage y_storage;
      typename VectorXsStorage::RefMapType y;

      /// \brief Temporary vector for computations.
      VectorXsStorage rhs_storage;
      typename VectorXsStorage::RefMapType rhs;

      /// \brief Temporary vector for computations.
      VectorXsStorage tmp_storage;
      typename VectorXsStorage::RefMapType tmp;
    }; // struct PGSSolverWorkspaceTpl
  } // namespace internal

} // namespace pinocchio

#include "pinocchio/algorithm/solvers/pgs-solver.hxx"

#endif // ifndef __pinocchio_algorithm_solvers_pgs_solver_hpp__
