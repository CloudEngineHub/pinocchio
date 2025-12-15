//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_pgs_solver_hpp__
#define __pinocchio_algorithm_pgs_solver_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/contact-solver-base.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include <boost/optional.hpp>
#include <limits>

namespace pinocchio
{

  template<typename Scalar>
  struct PGSConstraintSolverTpl;
  typedef PGSConstraintSolverTpl<context::Scalar> PGSConstraintSolver;

  ///
  /// \brief Settings for the PGS constraint solver loop.
  template<typename _Scalar>
  struct PGSSolverSettingsTpl : ConstraintSolverSettingsBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverSettingsBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;

    /// \brief Default constructor
    PGSSolverSettingsTpl(
      std::size_t max_iterations = 1000,
      Scalar tol_feasibility = Scalar(1e-6),
      Scalar tol_rel_feasibility = Scalar(1e-6),
      Scalar tol_complementarity = Scalar(1e-6),
      Scalar tol_rel_complementarity = Scalar(1e-6),
      bool solve_ncp = true,
      bool measure_timings = false,
      bool stat_record = false,
      Scalar over_relaxation = Scalar(1),
      std::optional<RefConstVectorXs> primal_guess = std::nullopt)
    : Base(
        max_iterations,
        tol_feasibility,
        tol_rel_feasibility,
        tol_complementarity,
        tol_rel_complementarity,
        solve_ncp,
        measure_timings,
        stat_record)
    , primal_guess(primal_guess)
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
    using Base::tol_feasibility;

    /// \brief Relative tolerance on the primal/dual feasibilityibility.
    using Base::tol_rel_feasibility;

    /// \brief Absolute tolerance of complementarity (duality complementarity).
    using Base::tol_complementarity;

    /// \brief Relative tolerance of complementarity (duality complementarity).
    using Base::tol_rel_complementarity;

    /// \brief Whether or not to solve the NCP. If set to solve, the equivalent CCP
    /// is solved.
    using Base::solve_ncp;

    /// \brief Measure solve timings
    using Base::measure_timings;

    /// \brief Record per iteration stats.
    using Base::stat_record;

    // ----------------------
    // Warmstart settings

    /// \brief Optional guess for the primal variable (impulses).
    std::optional<RefConstVectorXs> primal_guess;

    // ----------------------
    // PGS specific settings

    /// \brief Over-relaxation of PGS step. Default value is 1.
    Scalar over_relaxation;
  }; // struct PGSSolverSettingsTpl
  typedef PGSSolverSettingsTpl<context::Scalar> PGSSolverSettings;

  ///
  /// \brief Workspace for the PGS constraint solver.
  template<typename _Scalar>
  struct PGSSolverWorkspaceTpl
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    // TODO: use eigenstorage

    /// \brief Default constructor
    PGSSolverWorkspaceTpl()
    : problem_size(0)
    {
    }

    /// \brief Constructor given problem_size.
    PGSSolverWorkspaceTpl(std::size_t problem_size)
    : problem_size(problem_size)
    {
      resize(problem_size);
    }

    /// \brief Resize workspace vectors to problem size.
    void resize(std::size_t problem_size_)
    {
      problem_size = problem_size_;

      Eigen::Index np = static_cast<Eigen::Index>(problem_size);
      x.setZero(np);
      x_previous.setZero(np);
      y.setZero(np);
      tmp.setZero(np);
      rhs.setZero(np);
    }

    /// \brief Size of problem.
    std::size_t problem_size;

    /// \brief Primal variable (impulses) at current iteration.
    VectorXs x;

    /// \brief Primal variable (impulses) at previous iteration.
    VectorXs x_previous;

    /// \brief Dual variable (constraint velocities) at current iteration.
    VectorXs y;

    /// \brief Temporary vector for computations.
    VectorXs rhs;

    /// \brief Temporary vector for computations.
    VectorXs tmp;
  }; // struct PGSSolverWorkspaceTpl
  typedef PGSSolverWorkspaceTpl<context::Scalar> PGSSolverWorkspace;

  ///
  /// \brief Struct describing the solution of the PGS constraint solver
  /// after calling the `solve` method.
  template<typename _Scalar>
  struct PGSSolverSolutionTpl : ConstraintSolverSolutionBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverSolutionBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    using Base::nan;

    /// \brief Default constructor.
    PGSSolverSolutionTpl()
    : Base()
    , problem_size(0)
    {
    }

    /// \brief Reset the solution.
    void reset()
    {
      Base::reset();
      problem_size = 0;
    }

    /// \brief Resize the primal/dual vectors of the solution.
    void resize(std::size_t problem_size_)
    {
      problem_size = problem_size_;

      Eigen::Index np = static_cast<Eigen::Index>(problem_size);
      x.setZero(np);
      y.setZero(np);
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

    /// \brief Primal solution.
    VectorXs x;

    /// \brief Dual solution.
    VectorXs y;
  }; // struct PGSSolverSolutionTpl

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
  typedef PGSSolverStatsTpl<context::Scalar> PGSSolverStats;

  /// \brief Projected Gauss Siedel solver
  template<typename _Scalar>
  struct PGSConstraintSolverTpl : ConstraintSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;

    typedef PGSSolverWorkspaceTpl<Scalar> PGSSolverWorkspace;
    typedef PGSSolverSettingsTpl<Scalar> PGSSolverSettings;
    typedef PGSSolverSolutionTpl<Scalar> PGSSolverSolution;
    typedef PGSSolverStatsTpl<Scalar> PGSSolverStats;

    explicit PGSConstraintSolverTpl(std::size_t problem_size)
    : Base()
    , solution()
    , workspace(problem_size)
    , stats()
    {
    }

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_sets) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the contact problem.
    /// \param[in] g Free contact acceleration or velicity associted with the contact problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in] x Initial guess solution of the problem.
    /// \param[in] over_relax Optional over relaxation value, default to 1.
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
      const PGSSolverSettings & settings);

    ///
    /// \brief Solve the constrained problem composed of problem data (G,g,constraint_sets) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the contact problem.
    /// \param[in] g Free contact acceleration or velicity associted with the contact problem.
    /// \param[in] constraint_models Vector of constraint models.
    /// \param[in] x Initial guess solution of the problem.
    /// \param[in] over_relax Optional over relaxation value, default to 1.
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
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const PGSSolverSettings & settings)
    {
      return solve(delassus.derived().matrix(), g, constraint_models, constraint_datas, settings);
    }

#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL

    PGSSolverSolution solution;
    PGSSolverWorkspace workspace;
    PGSSolverStats stats;

  }; // struct PGSConstraintSolverTpl
} // namespace pinocchio

#include "pinocchio/algorithm/pgs-solver.hxx"

#endif // ifndef __pinocchio_algorithm_pgs_solver_hpp__
