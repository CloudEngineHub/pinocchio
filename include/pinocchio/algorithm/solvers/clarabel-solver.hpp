//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_solvers_clarabel_solver_hpp__
#define __pinocchio_algorithm_solvers_clarabel_solver_hpp__

#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT

  #include "pinocchio/algorithm/solvers/fwd.hpp"
  #include "pinocchio/algorithm/solvers/constraint-solver-base.hpp"

  #include "pinocchio/algorithm/delassus-operator-base.hpp"

  #include "pinocchio/algorithm/constraints/fwd.hpp"
  #include "pinocchio/math/fwd.hpp"

  #include "pinocchio/container/eigen-storage.hpp"

  #include <clarabel.hpp>

  #include <optional>
  #include <limits>

namespace pinocchio
{

  // Clarabel solver specific fwd declaration
  // see below for implementation
  namespace internal
  {
    template<typename Scalar, int Options>
    struct ClarabelSolverWorkspaceTpl;
  }

  ///
  /// \brief Clarabel-based cone solver wrapper for contact constraint resolution.
  ///
  /// This class provides a Pinocchio-compatible interface to the Clarabel
  /// interior point solver for solving contact problems with second-order cone
  /// constraints (friction cones).
  template<typename _Scalar, int _Options>
  struct ClarabelConstraintSolverTpl : ConstraintSolverBaseTpl<_Scalar>
  {

    using Scalar = _Scalar;
    static constexpr int Options = _Options;
    using Base = ConstraintSolverBaseTpl<_Scalar>;
    using VectorXs = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options>;
    using RefVectorXs = Eigen::Ref<VectorXs>;
    using RefConstVectorXs = Eigen::Ref<const VectorXs>;
    using ConstRefVectorXs = const Eigen::Ref<const VectorXs>;
    using MatrixXs = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor>;
    using Vector3s = Eigen::Matrix<Scalar, 3, 1, Options>;

    using ClarabelSolverResult = ClarabelSolverResultTpl<Scalar, Options>;
    using ClarabelSolverSettings = ClarabelSolverSettingsTpl<Scalar>;
    using ClarabelSolverWorkspace = internal::ClarabelSolverWorkspaceTpl<Scalar, Options>;

    /// \brief Default constructor
    explicit ClarabelConstraintSolverTpl(std::size_t problem_size = 0)
    : Base()
    , workspace_(problem_size)
    , is_valid_(false)
    {
    }

    ///
    /// \brief Solve the contact problem using Clarabel
    ///
    /// The problem is formulated as:
    ///   minimize    (1/2) lambda^T G lambda + g^T lambda
    ///   subject to  lambda in K
    ///
    /// where K is a product of second-order cones (friction cones) and other
    /// constraint sets.
    ///
    /// This method follows the same pattern as ADMM/PGS solvers in pinocchio.
    ///
    /// \param[in] delassus The Delassus operator G
    /// \param[in] g The drift term
    /// \param[in] constraint_models Vector of constraint models
    /// \param[in] constraint_datas Vector of constraint data
    /// \param[in] settings Settings for the Clarabel solver (max_iter,
    /// tolerances, etc.)
    /// \param[in/out] result Solution and warmstart data. Contains primal_guess,
    /// dual_guess, preconditioner for warmstarting.
    ///                       After solving, contains the solution in x, y, z,
    ///                       desaxce.
    ///
    /// \returns True if the problem has converged.
    ///
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
      const ClarabelSolverSettings & settings,
      ClarabelSolverResult & result);

    /// \brief Reset the constraint solver as if it has never run.
    void reset()
    {
      workspace_.reset();
      is_valid_ = false;
    }

    /// \brief Returns true if solver is in a valid state (it has solved a
    /// constraint problem).
    bool isValid() const
    {
      return is_valid_;
    }

  #ifdef PINOCCHIO_WITH_COLLISION
    /// \brief Timer for the `solve` method
    using Base::timer;
  #endif // PINOCCHIO_WITH_COLLISION

  protected:
    /// \brief Workspace of the Clarabel solver.
    /// This is an internal of the solver and is not meant to be accessed by
    /// users.
    ClarabelSolverWorkspace workspace_;

    /// \brief Flag to check whether or not the solver is in a reset state.
    bool is_valid_;
  };

  ///
  /// \brief Settings for the Clarabel constraint solver.
  ///
  /// This struct configures the Clarabel interior point solver for contact
  /// problems. Follows the same pattern as pinocchio::ADMMSolverSettingsTpl.
  ///
  template<typename _Scalar>
  struct ClarabelSolverSettingsTpl : ConstraintSolverSettingsBaseTpl<_Scalar>
  {
    using Scalar = _Scalar;
    using Base = ConstraintSolverSettingsBaseTpl<_Scalar>;

    /// \brief Default constructor.
    ClarabelSolverSettingsTpl(
      std::size_t max_iterations = 1000,
      Scalar absolute_feasibility_tol = Scalar(1e-8),
      Scalar relative_feasibility_tol = Scalar(1e-8),
      Scalar absolute_complementarity_tol = Scalar(1e-8),
      Scalar relative_complementarity_tol = Scalar(1e-8),
      bool solve_ncp = true,
      bool measure_timings = false,
      bool stat_record = false,
      std::size_t max_ncp_loops = 25,
      Scalar tol_ktratio = Scalar(1e-8),
      bool verbose = false)
    : Base(
        max_iterations,
        absolute_feasibility_tol,
        relative_feasibility_tol,
        absolute_complementarity_tol,
        relative_complementarity_tol,
        solve_ncp,
        measure_timings,
        stat_record)
    , max_ncp_loops(max_ncp_loops)
    , tol_ktratio(tol_ktratio)
    , verbose(verbose)
    {
    }

    /// \brief Throws if settings are not valid
    void checkValidity() const
    {
      Base::checkValidity();
      PINOCCHIO_CHECK_INPUT_ARGUMENT(tol_ktratio >= Scalar(0), "tol_ktratio should be >= 0.");
    }

    // ----------------------
    // General settings (inherited from base)

    /// \brief Maximum number of iterations of the solver.
    using Base::max_iterations;

    /// \brief Tolerance on the primal/dual feasibility.
    using Base::absolute_feasibility_tol;

    /// \brief Relative tolerance on the primal/dual feasibility.
    using Base::relative_feasibility_tol;

    /// \brief Absolute tolerance on the complementarity (duality gap).
    using Base::absolute_complementarity_tol;

    /// \brief Relative tolerance on the complementarity (duality gap).
    using Base::relative_complementarity_tol;

    /// \brief Whether or not to solve the NCP. If set to false, the equivalent
    /// CCP is solved.
    using Base::solve_ncp;

    /// \brief Measure solve timings
    using Base::measure_timings;

    // ----------------------
    // Clarabel specific settings

    /// \brief Maximum number of calls to the clarabel solver when solving the NCP problem.
    std::size_t max_ncp_loops;

    /// \brief Tolerance for KKT conditions (Clarabel-specific)
    Scalar tol_ktratio;

    /// \brief Whether to enable verbose output (Clarabel-specific)
    bool verbose;
  }; // struct ClarabelSolverSettingsTpl

  ///
  /// \brief Result container for Clarabel constraint solver.
  template<typename _Scalar, int _Options>
  struct ClarabelSolverResultTpl : ConstraintSolverResultBaseTpl<_Scalar>
  {
    using Scalar = _Scalar;
    static constexpr int Options = _Options;
    using Base = ConstraintSolverResultBaseTpl<_Scalar>;
    using VectorXs = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options>;
    using VectorXsStorage = EigenStorageTpl<VectorXs>;

    using Base::isValid;

    /// \brief Default constructor
    ClarabelSolverResultTpl()
    : Base()
    , problem_size(0)
    , x(x_storage.map())
    , z(z_storage.map())
    , desaxce(desaxce_storage.map())
    , sigma(sigma_storage.map())
    {
      reset(problem_size);
    }

    /// \brief Reset the results.
    void reset(std::size_t problem_size_ = 0)
    {
      Base::reset();
      problem_size = problem_size_;
      resize(problem_size);

      // Set solution to zero - solver has not run
      x.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      z.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      desaxce.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      sigma.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
    }

    /// \brief Resize the primal/dual/desaxce vectors.
    void resize(std::size_t problem_size_)
    {
      problem_size = problem_size_;
      Eigen::Index np = static_cast<Eigen::Index>(problem_size);

      x_storage.resize(np);
      z_storage.resize(np);
      desaxce_storage.resize(np);
      sigma_storage.resize(np);
    }

    /// \brief Retrieve primal solution (constraint impulses).
    template<typename VectorLike>
    void retrievePrimalSolution(const Eigen::MatrixBase<VectorLike> & primal_solution_) const
    {
      auto & primal_solution = primal_solution_.const_cast_derived();
      primal_solution = x;
    }

    /// \brief Retrieve dual solution (scaled constraint velocities)
    template<typename VectorLike>
    void retrieveDualSolution(const Eigen::MatrixBase<VectorLike> & dual_solution_) const
    {
      auto & dual_solution = dual_solution_.const_cast_derived();
      dual_solution = z;
    }

    /// \brief Retrieve De Saxcé term
    template<typename VectorLike>
    void retrieveDesaxceTerm(const Eigen::MatrixBase<VectorLike> & desaxce_term_) const
    {
      auto & desaxce_term = desaxce_term_.const_cast_derived();
      desaxce_term = desaxce;
    }

    /// \brief Retrieve constraints impulses.
    template<typename VectorLike>
    void
    retrieveConstraintImpulses(const Eigen::MatrixBase<VectorLike> & constraint_impulses_) const
    {
      auto & constraint_impulses = constraint_impulses_.const_cast_derived();
      constraint_impulses = x;
    }

    /// \brief Retrieve dual solution (scaled constraint velocities)
    template<typename VectorLike>
    void
    retrieveConstraintVelocities(const Eigen::MatrixBase<VectorLike> & constraint_velocities_) const
    {
      auto & constraint_velocities = constraint_velocities_.const_cast_derived();
      constraint_velocities = sigma;
    }

    // ----------------------
    // Inherited members from base class

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

    /// \brief Problem size
    std::size_t problem_size;

    // ----------------------
    // Solution - output of the solver

    /// \brief Primal solution vector (constraints impulses).
    VectorXsStorage x_storage;
    typename VectorXsStorage::RefMapType x;

    /// \brief Dual solution vector (scaled constraint velocities).
    VectorXsStorage z_storage;
    typename VectorXsStorage::RefMapType z;

    /// \brief De Saxcé shift term for NCP problems.
    VectorXsStorage desaxce_storage;
    typename VectorXsStorage::RefMapType desaxce;

    /// \brief Constraints velocities.
    VectorXsStorage sigma_storage;
    typename VectorXsStorage::RefMapType sigma;

  }; // struct ClarabelSolverResultTpl

  namespace internal
  {
    ///
    /// \brief Workspace for the Clarabel constraint solver.
    template<typename _Scalar, int _Options>
    struct ClarabelSolverWorkspaceTpl
    {
      using Scalar = _Scalar;
      static constexpr int Options = _Options;
      using VectorXs = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options>;
      using VectorXsStorage = EigenStorageTpl<VectorXs>;

      /// \brief Constructor given problem_size.
      explicit ClarabelSolverWorkspaceTpl(std::size_t problem_size = 0)
      : problem_size(problem_size)
      , clarabel_solver(nullptr)
      , sigma(sigma_storage.map())
      , desaxce(desaxce_storage.map())
      , prev_desaxce(prev_desaxce_storage.map())
      {
        resize(problem_size);
      }

      /// \brief Reset the workspace.
      void reset(std::size_t problem_size_ = 0)
      {
        problem_size = problem_size_;
        clarabel_solver.reset();
        clarabel_cones.clear();
        resize(problem_size);

  #ifndef NDEBUG
        // for debugging purposes
        sigma.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        desaxce.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
        prev_desaxce.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
  #endif
      }

      /// \brief Resize workspace vectors to problem size.
      void resize(std::size_t problem_size_)
      {
        problem_size = problem_size_;
        const Eigen::Index np = static_cast<Eigen::Index>(problem_size);

        sigma_storage.resize(np);
        desaxce_storage.resize(np);
        prev_desaxce_storage.resize(np);
      }

      /// \brief Size of problem.
      std::size_t problem_size;

      /// \brief Internal clarabel solver
      std::shared_ptr<::clarabel::DefaultSolver<Scalar>> clarabel_solver;

      /// \brief Internal cones for the clarabel solver
      std::vector<::clarabel::SupportedConeT<Scalar>> clarabel_cones;

      /// \brief Constraint velocities.
      VectorXsStorage sigma_storage;
      typename VectorXsStorage::RefMapType sigma;

      /// \brief De Saxcé shift.
      VectorXsStorage desaxce_storage;
      typename VectorXsStorage::RefMapType desaxce;

      /// \brief Previous NCP loop De Saxcé shift
      VectorXsStorage prev_desaxce_storage;
      typename VectorXsStorage::RefMapType prev_desaxce;
    }; // struct ClarabelSolverWorkspaceTpl

  } // namespace internal

} // namespace pinocchio

#endif // PINOCCHIO_WITH_CLARABEL_SUPPORT

// Include implementation
#include "pinocchio/src/algorithm/solvers/clarabel-solver.hxx"

#if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #include "pinocchio/algorithm/solvers/clarabel-solver.txx"
#endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // __pinocchio_algorithm_solvers_clarabel_solver_hpp__
