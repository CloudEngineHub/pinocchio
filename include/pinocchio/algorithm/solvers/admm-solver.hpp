//
// Copyright (c) 2022-2025 INRIA
//

#ifndef __pinocchio_algorithm_solvers_admm_solver_hpp__
#define __pinocchio_algorithm_solvers_admm_solver_hpp__

#include "pinocchio/algorithm/solvers/fwd.hpp"
#include "pinocchio/algorithm/solvers/constraint-solver-base.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"
#include "pinocchio/math/eigenvalues.hpp"

#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/solvers/anderson-acceleration.hpp"

#include "pinocchio/math/lanczos-decomposition.hpp"

#include "pinocchio/algorithm/diagonal-preconditioner.hpp"

#include "pinocchio/container/eigen-storage.hpp"

#include <optional>
#include <limits>

namespace pinocchio
{
  // fwd declarations for ADMM-internal structs.
  // user-api structs are fwd delclared in solvers/fwd.hpp.
  // see bottom of file for definitions
  namespace internal
  {
    template<typename Scalar>
    struct ADMMSolverWorkspaceTpl;

    template<typename Scalar>
    struct ADMMSpectralUpdateRuleTpl;

    template<typename Scalar>
    struct ADMMLinearUpdateRuleTpl;

    template<typename Scalar>
    struct ADMMOSQPUpdateRuleTpl;

    /// \brief Container of ADMM rho update rules.
    template<typename Scalar>
    union ADMMUpdateRuleContainerTpl {
      ADMMUpdateRuleContainerTpl()
      : dummy() {};
      ADMMSpectralUpdateRuleTpl<Scalar> spectral_rule;
      ADMMOSQPUpdateRuleTpl<Scalar> osqp_rule;
      ADMMLinearUpdateRuleTpl<Scalar> linear_rule;

    protected:
      struct Dummy
      {
        Dummy() {};
      };

      Dummy dummy{};
    }; // struct ADMMUpdateRuleContainerTpl
  } // namespace internal

  ///
  /// \brief ADMM rho update rule.
  /// SPECTRAL: if primal/dual ratio met, multiply rho by a spectral factor.
  /// OSQP: if primal/dual ratio met, multiply rho by sqrt(primal_feas/dual_feas).
  /// LINEAR: if primal/dual ratio met, multiply rho by a linear factor.
  /// CONSTANT: keep the same rho.
  enum class ADMMUpdateRule : char
  {
    SPECTRAL = 'S',
    OSQP = 'O',
    LINEAR = 'L',
    CONSTANT = 'C',
  };

  ///
  /// \brief ADMM proximal policy.
  /// MANUAL: mu_prox is constant and manually set. It is scaled by tau_prox
  /// AUTOMATIC: mu_prox is always set to rho.
  ///
  /// \note mu_prox is always scaled by tau_prox.
  enum class ADMMProximalRule : char
  {
    MANUAL = 'M',
    AUTOMATIC = 'A',
  };

  ///
  /// \brief ADMM constraint solver.
  ///
  /// The solver solves the following CPP/NCP problem:
  /// `min_x x^T G x + g s.t. x \in C`,
  /// where `G` is the delassus matrix, `g` is the constraint velocities without any constraint
  /// forces acting on the system and `C` are the constraint sets. If the `g` term is augmented with
  /// the DeSaxce term, the problem becomes an NCP.
  template<typename _Scalar>
  struct PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
    ADMMConstraintSolverTpl : ConstraintSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverBaseTpl<_Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<VectorXs> RefVectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;
    typedef const Eigen::Ref<const VectorXs> ConstRefVectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    typedef ADMMSolverResultTpl<Scalar> ADMMSolverResult;
    typedef internal::ADMMSolverWorkspaceTpl<Scalar> ADMMSolverWorkspace;
    typedef ADMMSolverSettingsTpl<Scalar> ADMMSolverSettings;
    typedef ADMMSolverStatsTpl<Scalar> ADMMSolverStats;

    typedef internal::ADMMSpectralUpdateRuleTpl<Scalar> ADMMSpectralUpdateRule;
    typedef internal::ADMMOSQPUpdateRuleTpl<Scalar> ADMMOSQPUpdateRule;
    typedef internal::ADMMLinearUpdateRuleTpl<Scalar> ADMMLinearUpdateRule;
    typedef internal::ADMMUpdateRuleContainerTpl<Scalar> ADMMUpdateRuleContainer;

    /// \brief Default constructor.
    explicit ADMMConstraintSolverTpl(std::size_t problem_size = 0)
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
    /// \param[in] settings Settings for the ADMM solver.
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
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const ADMMSolverSettings & settings,
      ADMMSolverResult & result);

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

    /// \brief Per-iteration stats of the ADMM solver.
    ADMMSolverStats stats;

  protected:
    /// \brief Workspace of the ADMM solver.
    /// This is an internal of the solver and is not meant to be accessed by
    /// users.
    ADMMSolverWorkspace m_workspace;

    /// \brief Flag to check whether or not the solver is in a reset state.
    /// If not, its stats are valid.
    bool m_is_valid;

    /// \brief Compute largest eigen value of delassus.
    template<typename DelassusDerived>
    static Scalar computeDelassusLargestEigenvalue(
      const DelassusOperatorBase<DelassusDerived> & delassus, ADMMSolverWorkspace & workspace);

    /// \brief Retrieve primal and/or dual guesses from settings and result's warmstarts.
    template<
      typename DelassusDerived,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    static void retrievePrimalDualGuess(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const ADMMSolverSettings & settings,
      const ADMMSolverResult & result,
      ADMMSolverWorkspace & workspace);

    /// \brief Retrieve rho parameters guesses from settings and result's warmstarts.
    template<typename DelassusDerived>
    static void retrieveRhoGuess(
      const DelassusOperatorBase<DelassusDerived> & delassus,
      const ADMMSolverSettings & settings,
      const ADMMSolverResult & result,
      ADMMSolverWorkspace & workspace);
  }; // struct ADMMConstraintSolverTpl

  ///
  /// \brief Settings for the ADMM constraint solver loop.
  template<typename _Scalar>
  struct ADMMSolverSettingsTpl : ConstraintSolverSettingsBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverSettingsBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;

    /// \brief Default constructor
    ADMMSolverSettingsTpl(
      std::size_t max_iterations = 1000,
      Scalar absolute_feasibility_tol = Scalar(1e-6),
      Scalar relative_feasibility_tol = Scalar(1e-6),
      Scalar absolute_complementarity_tol = Scalar(1e-6),
      Scalar relative_complementarity_tol = Scalar(1e-6),
      bool solve_ncp = true,
      bool measure_timings = false,
      bool stat_record = false,
      std::optional<Scalar> rho_init = std::nullopt,
      bool warmstart_rho_with_previous_result = false,
      ADMMUpdateRule admm_update_rule = ADMMUpdateRule::OSQP,
      ADMMProximalRule admm_proximal_rule = ADMMProximalRule::MANUAL,
      Scalar mu_prox = Scalar(1e-6),
      Scalar tau_prox = Scalar(1),
      Scalar tau = Scalar(0.5),
      Scalar ratio_primal_dual = Scalar(10),
      Scalar dual_momentum = Scalar(0),
      Scalar rho_update_ratio = Scalar(0),
      std::size_t rho_min_update_frequency = 1,
      Scalar rho_momentum = Scalar(0),
      Scalar rho_min = Scalar(1e-6),
      Scalar rho_max = Scalar(1e6),
      Scalar spectral_rho_power_init = Scalar(0.2),
      Scalar spectral_rho_power_factor = Scalar(0.05),
      Scalar linear_update_rule_factor = Scalar(2),
      std::size_t lanczos_size = std::numeric_limits<int>::max(),
      std::size_t max_delassus_decomposition_updates = std::numeric_limits<int>::max(),
      std::size_t anderson_capacity = 0)
    : Base(
        max_iterations,
        absolute_feasibility_tol,
        relative_feasibility_tol,
        absolute_complementarity_tol,
        relative_complementarity_tol,
        solve_ncp,
        measure_timings,
        stat_record)
    , rho_init(rho_init)
    , warmstart_rho_with_previous_result(warmstart_rho_with_previous_result)
    , admm_update_rule(admm_update_rule)
    , admm_proximal_rule(admm_proximal_rule)
    , mu_prox(mu_prox)
    , tau_prox(tau_prox)
    , tau(tau)
    , ratio_primal_dual(ratio_primal_dual)
    , dual_momentum(dual_momentum)
    , rho_update_ratio(rho_update_ratio)
    , rho_min_update_frequency(rho_min_update_frequency)
    , rho_momentum(rho_momentum)
    , rho_min(rho_min)
    , rho_max(rho_max)
    , spectral_rho_power_init(spectral_rho_power_init)
    , spectral_rho_power_factor(spectral_rho_power_factor)
    , linear_update_rule_factor(linear_update_rule_factor)
    , lanczos_size(lanczos_size)
    , max_delassus_decomposition_updates(max_delassus_decomposition_updates)
    , anderson_capacity(anderson_capacity)
    {
    }

    /// \brief Throws if settings are not valid.
    void checkValidity() const
    {
      Base::checkValidity();
      if (rho_init)
      {
        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          rho_init.value() > Scalar(0), "rho_init should be none or > 0.");
      }
      PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_prox > Scalar(0), "mu_prox should be > 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        tau_prox <= Scalar(1) && tau_prox > Scalar(0), "tau_prox should lie in ]0,1].");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        tau <= Scalar(1) && tau > Scalar(0), "tau should lie in ]0,1].");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        ratio_primal_dual > Scalar(0), "ratio_primal_dual should be > 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(dual_momentum >= Scalar(0), "dual_momentum should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        rho_momentum >= Scalar(0) && rho_momentum <= Scalar(1),
        "rho_momentum should be in [0, 1].");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        linear_update_rule_factor >= Scalar(0), "linear_update_rule_factor should be >= 0.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        rho_min >= Scalar(0) && rho_max >= Scalar(0) && rho_min <= rho_max,
        "rho_min and rho_max don't verify 0 <= rho_min <= rho_max");
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

    /// \brief Whether or not to solve the NCP. If set to false, the equivalent CCP
    /// is solved. Default is true.
    using Base::solve_ncp;

    /// \brief Measure solve timings
    using Base::measure_timings;

    /// \brief Record per iteration stats.
    using Base::stat_record;

    // ----------------------
    // ADMM specific settings

    /// \brief Initial value of rho parameter.
    /// If set to boost::none, the initial rho will be computed by estimating
    /// the largest eigenvalue of the Delassus.
    std::optional<Scalar> rho_init;

    /// \brief Whether or not to warmstart rho with previous result.
    /// If set to true, the `rho_init` will be bypassed by the value of rho
    /// stored in the solver's solution.
    /// \note This setting is only effective if a previous solution is valid.
    bool warmstart_rho_with_previous_result;

    /// \brief Update rule for the rho admm term.
    ADMMUpdateRule admm_update_rule;

    /// \brief Update rule for the primal proximal term.
    ADMMProximalRule admm_proximal_rule;

    /// \brief Value of the proximal term when `admm_proximal_rule` is `MANUAL`.
    /// When `admm_proximal_rule` is `AUTOMATIC`, the proximal term follows the rho term.
    Scalar mu_prox;

    /// \brief Scaling factor in front of the primal proximal term.
    Scalar tau_prox;

    /// \brief Scaling factor in front of the rho ADMM term.
    Scalar tau;

    /// \brief Value of the primal/dual ratio beyond/below which a rho update is considered.
    /// If the primal/dual ratio is not big/low enough, rho is kept constant.
    Scalar ratio_primal_dual;

    /// \brief Momentum on the dual variable. 0 means no momentum.
    Scalar dual_momentum;

    /// \brief Ratio w.r.t previous rho value beyond/below which rho is updated.
    /// In essence, if rho has not changed enough compared to the previous rho value,
    /// its value stays unchanged. This is a knob to prevent updating rho too much.
    Scalar rho_update_ratio;

    /// \brief How many iterations before rho can be updated again. 1 means rho can be updated.
    /// every iterations (it does not mean it is updated every iterations).
    std::size_t rho_min_update_frequency;

    /// \brief Momentum on rho value
    Scalar rho_momentum;

    /// \brief Minimum value that rho can take.
    /// During 'solve', rho is clamped between `rho_min` and `rho_max`.
    Scalar rho_min;

    /// \brief Maximum value that rho can take.
    /// During 'solve', rho is clamped between `rho_min` and `rho_max`.
    Scalar rho_max;

    /// \brief Initial value of the rho power in the SPECTRAL update rule.
    Scalar spectral_rho_power_init;

    /// \brief Power factor used to update rho in the SPECTRAL update rule.
    Scalar spectral_rho_power_factor;

    /// \brief Value by which rho is multiplied/divided in the LINEAR update rule.
    Scalar linear_update_rule_factor;

    /// \brief Size of the lanczos decomposition.
    /// Higher values lead to more precise estimation of the Delassus' maximum eigenvalue.
    std::size_t lanczos_size;

    /// \brief Maximum number of delassus decomposition updates.
    /// Once this number is reached, rho and mu_prox are kept constant.
    std::size_t max_delassus_decomposition_updates;

    /// \brief Size of the andersion history used to fit the anderson linear system.
    /// If set to 0 or 1, no anderson acceleration is used.
    std::size_t anderson_capacity;

  }; // struct ADMMSolverSettingsTpl

  ///
  /// \brief Struct describing the solution of the ADMM constraint solver
  /// after calling the `solve` method.
  /// Also contains the warmstart of the solution to the constraint problem.
  template<typename _Scalar>
  struct ADMMSolverResultTpl : ConstraintSolverResultBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverResultBaseTpl<Scalar> Base;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;
    typedef EigenStorageTpl<VectorXs> VectorXsStorage;

    // make ADMM solver a friend so that it can use `makeValid`
    friend struct ADMMConstraintSolverTpl<Scalar>;

    using Base::isValid;
    /// \brief Default constructor.
    ADMMSolverResultTpl()
    : Base()
    , problem_size(0)
    , delassus_decomposition_update_count(0)
    , preconditioner(std::nullopt)
    , primal_guess(std::nullopt)
    , dual_guess(std::nullopt)
    , rho(std::numeric_limits<Scalar>::quiet_NaN())
    , spectral_rho_power(std::numeric_limits<Scalar>::quiet_NaN())
    , mu_prox(std::numeric_limits<Scalar>::quiet_NaN())
    , x(x_storage.map())
    , y(y_storage.map())
    , z(z_storage.map())
    , desaxce(desaxce_storage.map())
    {
    }

    /// \brief Reset the results.
    void reset(std::size_t problem_size_ = 0)
    {
      Base::reset();
      problem_size = problem_size_;
      delassus_decomposition_update_count = 0;

      preconditioner.reset();
      primal_guess.reset();
      dual_guess.reset();

      rho = std::numeric_limits<Scalar>::quiet_NaN();
      spectral_rho_power = std::numeric_limits<Scalar>::quiet_NaN();
      mu_prox = std::numeric_limits<Scalar>::quiet_NaN();

      resize(problem_size);

      // set solution to nan - solver has not run
      x.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      y.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      z.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
      desaxce.setConstant(std::numeric_limits<Scalar>::quiet_NaN());
    }

    /// \brief Resize the primal/dual/desaxce vectors of the solution.
    void resize(std::size_t problem_size_)
    {
      problem_size = problem_size_;

      Eigen::Index np = static_cast<Eigen::Index>(problem_size);
      x_storage.resize(np);
      y_storage.resize(np);
      z_storage.resize(np);
      desaxce_storage.resize(np);
    }

    /// \brief Retrieve primal solution.
    template<typename VectorLike>
    void retrievePrimalSolution(const Eigen::MatrixBase<VectorLike> & primal_solution_) const
    {
      auto & primal_solution = primal_solution_.const_cast_derived();
      primal_solution = y;
    }

    /// \brief Retrieve dual solution.
    template<typename VectorLike>
    void retrieveDualSolution(const Eigen::MatrixBase<VectorLike> & dual_solution_) const
    {
      auto & dual_solution = dual_solution_.const_cast_derived();
      dual_solution = z;
    }

    /// \brief Retrieve Desaxce term of solution
    template<typename VectorLike>
    void retrieveDesaxceTerm(const Eigen::MatrixBase<VectorLike> & desaxce_term_) const
    {
      auto & desaxce_term = desaxce_term_.const_cast_derived();
      desaxce_term = desaxce;
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

    /// \brief Number of delassus decompositions.
    std::size_t delassus_decomposition_update_count;

    // ----------------------
    // Solution warmstart

    /// \brief Optional preconditioner. Not used yet.
    std::optional<RefConstVectorXs> preconditioner;

    /// \brief Optional guess for the primal variable (impulses).
    std::optional<RefConstVectorXs> primal_guess;

    /// \brief Optional guess for the dual variable (velocities).
    std::optional<RefConstVectorXs> dual_guess;

    // ----------------------
    // Solution - output of the solver

    /// \brief Value of ADMM rho term.
    Scalar rho;

    /// \brief Value of ADMM spectral rule rho power at the solution.
    /// This is relevant only if `SPECTRAL` was selected as an update rule.
    Scalar spectral_rho_power;

    /// \brief Value of ADMM proximal term.
    Scalar mu_prox;

    /// \brief Non-projected primal solution.
    /// \note Order of storage/map declaration is important!
    /// First declare the storage, then the map, otherwise map will point to nothing.
    VectorXsStorage x_storage;
    typename VectorXsStorage::RefMapType x;

    /// \brief Primal solution projected onto constraints.
    VectorXsStorage y_storage;
    typename VectorXsStorage::RefMapType y;

    /// \brief Dual solution.
    VectorXsStorage z_storage;
    typename VectorXsStorage::RefMapType z;

    /// \brief Desaxce term of the solution
    VectorXsStorage desaxce_storage;
    typename VectorXsStorage::RefMapType desaxce;
  }; // struct ADMMSolverResultTpl

  ///
  /// \brief Struct to track per iteration progress of ADMM constraint solver.
  template<typename _Scalar>
  struct ADMMSolverStatsTpl : ConstraintSolverStatsBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverStatsBaseTpl<Scalar> Base;

    /// \brief Default constructor.
    ADMMSolverStatsTpl()
    : Base()
    , delassus_decomposition_update_count(0)
    {
    }

    /// \brief Constructor given a maximum iteration of the solver.
    explicit ADMMSolverStatsTpl(std::size_t max_iterations)
    : Base(max_iterations)
    , delassus_decomposition_update_count(0)
    {
      reserve(max_iterations);
    }

    /// \brief Reserve enough storage for max_it iterations.
    void reserve(std::size_t max_iterations)
    {
      rho.reserve(size_t(max_iterations));
      mu_prox.reserve(size_t(max_iterations));
      anderson_size.reserve(size_t(max_iterations));
      linear_system_residual.reserve(size_t(max_iterations));
      linear_system_consistency.reserve(size_t(max_iterations));
    }

    /// \brief Reset stats.
    void reset()
    {
      Base::reset();
      rho.clear();
      mu_prox.clear();
      anderson_size.clear();
      linear_system_residual.clear();
      linear_system_consistency.clear();
      delassus_decomposition_update_count = 0;
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

    /// \brief History of rho values.
    std::vector<Scalar> rho;

    /// \brief History of mu_prox values.
    std::vector<Scalar> mu_prox;

    /// \brief History of Anderson size.
    std::vector<std::size_t> anderson_size;

    /// \brief History of linear system residual.
    std::vector<Scalar> linear_system_residual;

    /// \brief History of linear system consistency
    std::vector<Scalar> linear_system_consistency;

    /// \brief Number of Delassus decomposition updates.
    std::size_t delassus_decomposition_update_count;
  }; // struct ADMMSolverStatsTpl

  namespace internal
  {
    ///
    /// \brief Workspace for the ADMM constraint solver.
    template<typename _Scalar>
    struct ADMMSolverWorkspaceTpl
    {
      typedef _Scalar Scalar;
      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
      typedef EigenStorageTpl<VectorXs> VectorXsStorage;
      typedef LanczosDecompositionTpl<MatrixXs> LanczosDecomposition;
      typedef AndersonAccelerationTpl<Scalar> AndersonAcceleration;

      /// \brief Constructor given problem_size, lanczos_size and anderson_capacity.
      ADMMSolverWorkspaceTpl(
        std::size_t problem_size = 0, //
        std::size_t lanczos_size = 2, //
        std::size_t anderson_capacity = 0)
      : problem_size(problem_size)
      , lanczos_size(lanczos_size)
      , anderson_capacity(anderson_capacity)
      , delassus_decomposition_update_count(0)
      , delassus_smallest_eigenvalue(std::nullopt)
      , delassus_largest_eigenvalue(std::nullopt)
      , rho(std::numeric_limits<Scalar>::quiet_NaN())
      , spectral_rho_power(std::numeric_limits<Scalar>::quiet_NaN())
      , mu_prox(std::numeric_limits<Scalar>::quiet_NaN())
      , lanczos_decomposition(
          static_cast<Eigen::Index>(math::max(std::size_t(2), problem_size)), //
          static_cast<Eigen::Index>(
            math::max(std::size_t(2), math::min(problem_size, lanczos_size))))
      , anderson_history(problem_size, anderson_capacity)
      , x(x_storage.map())
      , x_previous(x_previous_storage.map())
      , x_anderson(x_anderson_storage.map())
      , y(y_storage.map())
      , y_previous(y_previous_storage.map())
      , z(z_storage.map())
      , z_previous(z_previous_storage.map())
      , z_anderson(z_anderson_storage.map())
      , desaxce(desaxce_storage.map())
      , rhs(rhs_storage.map())
      , tmp(tmp_storage.map())
      , primal_feasibility_vector(primal_feasibility_vector_storage.map())
      , anderson_primal_feasibility_vector(anderson_primal_feasibility_vector_storage.map())
      , dual_feasibility_vector(dual_feasibility_vector_storage.map())
      {
        resize(problem_size, lanczos_size, anderson_capacity);
      }

      /// \brief Reset the workspace.
      void reset(
        std::size_t problem_size_ = 0, //
        std::size_t lanczos_size_ = 2, //
        std::size_t anderson_capacity_ = 0)
      {
        problem_size = problem_size_;
        lanczos_size = math::max(std::size_t(2), math::min(problem_size, lanczos_size_));
        anderson_capacity = anderson_capacity_;

        delassus_decomposition_update_count = 0;
        delassus_smallest_eigenvalue.reset();
        delassus_largest_eigenvalue.reset();
        constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();

        rho = nan;
        spectral_rho_power = nan;
        mu_prox = nan;

        resize(problem_size, lanczos_size, anderson_capacity);

#ifndef NDEBUG
        // for debugging purposes
        x.setConstant(nan);
        x_previous.setConstant(nan);
        x_anderson.setConstant(nan);
        y.setConstant(nan);
        y_previous.setConstant(nan);
        z.setConstant(nan);
        z_previous.setConstant(nan);
        z_anderson.setConstant(nan);
        desaxce.setConstant(nan);
        rhs.setConstant(nan);
        tmp.setConstant(nan);
        primal_feasibility_vector.setConstant(nan);
        anderson_primal_feasibility_vector.setConstant(nan);
        dual_feasibility_vector.setConstant(nan);
#endif
      }

      /// \brief Resize workspace vectors and operators to problem sizes.
      void resize(
        std::size_t problem_size_, //
        std::size_t lanczos_size_, //
        std::size_t anderson_capacity_)
      {
        problem_size = problem_size_;
        lanczos_size = math::max(std::size_t(2), math::min(problem_size, lanczos_size_));
        anderson_capacity = anderson_capacity_;

        const Eigen::Index np = static_cast<Eigen::Index>(problem_size);
        x_storage.resize(np);
        x_anderson_storage.resize(np);
        y_storage.resize(np);
        x_previous_storage.resize(np);
        y_previous_storage.resize(np);
        z_storage.resize(np);
        z_anderson_storage.resize(np);
        z_previous_storage.resize(np);
        desaxce_storage.resize(np);
        rhs_storage.resize(np);
        tmp_storage.resize(np);
        primal_feasibility_vector_storage.resize(np);
        anderson_primal_feasibility_vector_storage.resize(np);
        dual_feasibility_vector_storage.resize(np);

        // resize lanczos
        const std::size_t lanczos_problem_size = math::max(std::size_t(2), problem_size);
        if (
          lanczos_decomposition.size() != static_cast<Eigen::Index>(lanczos_problem_size)
          || lanczos_decomposition.decompositionSize() != static_cast<Eigen::Index>(lanczos_size))
        {
          lanczos_decomposition = LanczosDecomposition(
            static_cast<Eigen::Index>(lanczos_problem_size),
            static_cast<Eigen::Index>(lanczos_size));
        }

        // resize anderson
        anderson_history.reserve(problem_size, anderson_capacity);
      }

      /// \brief Size of problem.
      std::size_t problem_size;

      /// \brief Size of lanczos decomposition.
      std::size_t lanczos_size;

      /// \brief Size of anderson history capacity.
      std::size_t anderson_capacity;

      /// \brief Number of delassus decompositions.
      std::size_t delassus_decomposition_update_count;

      /// \brief Estimate of the smallest eigenvalue of the Delassus matrix.
      /// If its value is `std::nullopt`, then it has not been computed yet.
      std::optional<Scalar> delassus_smallest_eigenvalue;

      /// \brief Estimate of the largest eigenvalue of the Delassus matrix.
      /// If its value is `std::nullopt`, then it has not been computed yet.
      std::optional<Scalar> delassus_largest_eigenvalue;

      /// \brief Value of ADMM rho term at current iteration.
      Scalar rho;

      /// \brief Value of spectral rule rho power at current iteration.
      Scalar spectral_rho_power;

      /// \brief Value of primal proximal term at current iteration.
      Scalar mu_prox;

      /// \brief Lanczos decomposition used to estimate the maximum eigenvalue
      /// of the delassus. This allows to estimate the initial rho value if none
      /// is given as input to `solve`.
      /// The max eigenvalue is also needed in the SPECTRAL update rule.
      LanczosDecomposition lanczos_decomposition;

      /// \brief Anderson history used to fit the Anderson system.
      AndersonAcceleration anderson_history;

      /// \brief Primal variable at current iteration.
      VectorXsStorage x_storage;
      typename VectorXsStorage::RefMapType x;

      /// \brief Primal variable at previous iteration.
      VectorXsStorage x_previous_storage;
      typename VectorXsStorage::RefMapType x_previous;

      /// \brief Anderson accelerated primal variable at current iteration.
      VectorXsStorage x_anderson_storage;
      typename VectorXsStorage::RefMapType x_anderson;

      /// \brief Projected primal variable at current iteration.
      VectorXsStorage y_storage;
      typename VectorXsStorage::RefMapType y;

      /// \brief Projected primal variable at previous iteration.
      VectorXsStorage y_previous_storage;
      typename VectorXsStorage::RefMapType y_previous;

      /// \brief Dual variable at current iteration.
      VectorXsStorage z_storage;
      typename VectorXsStorage::RefMapType z;

      /// \brief Dual variable at previous iteration.
      VectorXsStorage z_previous_storage;
      typename VectorXsStorage::RefMapType z_previous;

      /// \brief Anderson accelerated dual variable at current iteration.
      VectorXsStorage z_anderson_storage;
      typename VectorXsStorage::RefMapType z_anderson;

      /// \brief Desaxce correction (always 0 when solving CCP) at current iteration.
      VectorXsStorage desaxce_storage;
      typename VectorXsStorage::RefMapType desaxce;

      /// \brief Temporary variable for solving linear system.
      VectorXsStorage rhs_storage;
      typename VectorXsStorage::RefMapType rhs;

      /// \brief Temporary variable for holding various vectors.
      VectorXsStorage tmp_storage;
      typename VectorXsStorage::RefMapType tmp;

      /// \brief Primal feasibility vector at current iteration.
      VectorXsStorage primal_feasibility_vector_storage;
      typename VectorXsStorage::RefMapType primal_feasibility_vector;

      /// \brief Anderson related primal feasibility vector at current iteration.
      VectorXsStorage anderson_primal_feasibility_vector_storage;
      typename VectorXsStorage::RefMapType anderson_primal_feasibility_vector;

      /// \brief Dual feasibility vector at current iteration.
      VectorXsStorage dual_feasibility_vector_storage;
      typename VectorXsStorage::RefMapType dual_feasibility_vector;
    }; // struct ADMMSolverWorkspaceTpl

    typedef ADMMSolverWorkspaceTpl<context::Scalar> ADMMSolverWorkspace;

    ///
    /// \brief Implementation of SPECTRAL ADMM update rule.
    template<typename _Scalar>
    struct ADMMSpectralUpdateRuleTpl
    {
      typedef _Scalar Scalar;

      /// \brief Default constructor.
      ADMMSpectralUpdateRuleTpl(
        const Scalar ratio_primal_dual,
        const Scalar L,
        const Scalar m,
        const Scalar rho_power_factor)
      : ratio_primal_dual(ratio_primal_dual)
      , rho_increment(std::pow(L / m, rho_power_factor))
      {
        PINOCCHIO_CHECK_INPUT_ARGUMENT(m > Scalar(0), "m should be positive.");
        PINOCCHIO_CHECK_INPUT_ARGUMENT(L > m, "L should be > m");
      }

      /// \brief Compute rho based on primal and dual feasibility.
      void eval(const Scalar primal_feasibility, const Scalar dual_feasibility, Scalar & rho) const
      {
        if (primal_feasibility > ratio_primal_dual * dual_feasibility)
        {
          rho *= rho_increment;
          //        rho *= math::pow(cond,rho_power_factor);
          //        rho_power += rho_power_factor;
        }
        else if (dual_feasibility > ratio_primal_dual * primal_feasibility)
        {
          rho /= rho_increment;
          //        rho *= math::pow(cond,-rho_power_factor);
          //        rho_power -= rho_power_factor;
        }
      }

      /// \brief Compute the penalty ADMM value from the current largest and lowest eigenvalues and
      /// the scaling spectral factor.
      static inline Scalar computeRho(const Scalar L, const Scalar m, const Scalar rho_power)
      {
        const Scalar cond = L / m;
        const Scalar rho = math::sqrt(L * m) * math::pow(cond, rho_power);
        return rho;
      }

      /// \brief Compute the scaling spectral factor of the ADMM penalty term from the current
      /// largest and lowest eigenvalues and the ADMM penalty term.
      static inline Scalar computeRhoPower(const Scalar L, const Scalar m, const Scalar rho)
      {
        const Scalar cond = L / m;
        const Scalar sqrt_L_m = math::sqrt(L * m);
        const Scalar rho_power = math::log(rho / sqrt_L_m) / math::log(cond);
        return rho_power;
      }

      /// \brief Primal/dual ratio below/above which rho is updated.
      Scalar ratio_primal_dual;

      /// \brief Factor by which rho gets divided/multiplied when primal/dual is below/above
      /// ratio_primal_dual.
      Scalar rho_increment;
    };

    ///
    /// \brief Implementation of LINEAR ADMM update rule.
    template<typename _Scalar>
    struct ADMMLinearUpdateRuleTpl
    {
      typedef _Scalar Scalar;

      /// \brief Constructor from primal/dual ratio and increase/decrease factors.
      ADMMLinearUpdateRuleTpl(
        const Scalar ratio_primal_dual, const Scalar increase_factor, const Scalar decrease_factor)
      : ratio_primal_dual(ratio_primal_dual)
      , increase_factor(increase_factor)
      , decrease_factor(decrease_factor)
      {
        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          increase_factor > Scalar(1), "increase_factor should be greater than one.");
        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          decrease_factor > Scalar(1), "decrease_factor should be greater than one.");
      }

      /// \brief Constructor from primal/dual ratio and factor (same one used to increase and
      /// decrease rho).
      ADMMLinearUpdateRuleTpl(const Scalar ratio_primal_dual, const Scalar factor)
      : ratio_primal_dual(ratio_primal_dual)
      , increase_factor(factor)
      , decrease_factor(factor)
      {
        PINOCCHIO_CHECK_INPUT_ARGUMENT(factor > Scalar(1), "factor should be greater than one.");
      }

      /// \brief Compute rho based on primal and dual feasibility.
      void eval(const Scalar primal_feasibility, const Scalar dual_feasibility, Scalar & rho) const
      {
        if (primal_feasibility > ratio_primal_dual * dual_feasibility)
        {
          rho *= increase_factor;
        }
        else if (dual_feasibility > ratio_primal_dual * primal_feasibility)
        {
          rho /= decrease_factor;
        }
      }

      /// \brief Primal/dual ratio below/above which rho is updated.
      Scalar ratio_primal_dual;

      /// \brief Factor by which rho gets multipled when primal/dual is above ratio.
      Scalar increase_factor;

      /// \brief Factor by which rho gets divided when primal/dual is below ratio.
      Scalar decrease_factor;
    };

    ///
    /// \brief Implementation of OSQP ADMM update rule.
    template<typename _Scalar>
    struct ADMMOSQPUpdateRuleTpl
    {
      typedef _Scalar Scalar;

      /// \brief Default constructor given a primal/dual ratio and a regularisation factor.
      ADMMOSQPUpdateRuleTpl(const Scalar ratio_primal_dual, const Scalar eps_reg = 1e-8)
      : ratio_primal_dual(ratio_primal_dual)
      , eps_reg(eps_reg)
      {
        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          ratio_primal_dual > Scalar(0), "ratio_primal_dual should be positive.");
        PINOCCHIO_CHECK_INPUT_ARGUMENT(eps_reg > Scalar(0), "eps_reg should be positive.");
      }

      /// \brief Compute rho based on primal and dual feasibility.
      void eval(const Scalar primal_feasibility, const Scalar dual_feasibility, Scalar & rho) const
      {
        if (
          primal_feasibility > ratio_primal_dual * dual_feasibility //
          || dual_feasibility > ratio_primal_dual * primal_feasibility)
        {
          rho *= std::sqrt(primal_feasibility / (dual_feasibility + eps_reg));
        }
      }

      /// \brief Primal/dual ratio below/above which rho is updated.
      Scalar ratio_primal_dual;

      /// \brief Regularisation factor to avoid dividing by 0.
      /// Default is 1e-8.
      Scalar eps_reg;
    };
  } // namespace internal

} // namespace pinocchio

#include "pinocchio/algorithm/solvers/admm-solver.hxx"

#endif // ifndef __pinocchio_algorithm_solvers_admm_solver_hpp__
