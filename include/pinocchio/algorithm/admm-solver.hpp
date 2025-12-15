//
// Copyright (c) 2022-2025 INRIA
//

#ifndef __pinocchio_algorithm_admm_solver_hpp__
#define __pinocchio_algorithm_admm_solver_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/comparison-operators.hpp"
#include "pinocchio/math/eigenvalues.hpp"

#include "pinocchio/algorithm/contact-solver-base.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/solvers/anderson-acceleration.hpp"

#include "pinocchio/math/lanczos-decomposition.hpp"

#include "pinocchio/algorithm/diagonal-preconditioner.hpp"

#include <optional>
#include <limits>

namespace pinocchio
{
  template<typename Scalar>
  struct ADMMConstraintSolverTpl;
  typedef ADMMConstraintSolverTpl<context::Scalar> ADMMConstraintSolver;

  ///
  /// \brief Implementation of SPECTRAL ADMM update rule.
  template<typename _Scalar>
  struct ADMMSpectralUpdateRuleTpl
  {
    typedef _Scalar Scalar;

    ADMMSpectralUpdateRuleTpl(
      const Scalar ratio_primal_dual, const Scalar L, const Scalar m, const Scalar rho_power_factor)
    : ratio_primal_dual(ratio_primal_dual)
    , rho_increment(std::pow(L / m, rho_power_factor))
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(m > Scalar(0), "m should be positive.");
    }

    Scalar getRatioPrimalDual() const
    {
      return ratio_primal_dual;
    }
    void setRatioPrimalDual(const Scalar ratio_primal_dual)
    {
      this->ratio_primal_dual = ratio_primal_dual;
    }

    Scalar getRhoIncrement() const
    {
      return rho_increment;
    }
    void setRhoIncrement(const Scalar cond, const Scalar rho_power_factor)
    {
      rho_increment = std::pow(cond, rho_power_factor);
    }

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

    /// \brief Compute the  scaling spectral factor of the ADMM penalty term from the current
    /// largest and lowest eigenvalues and the ADMM penalty term.
    static inline Scalar computeRhoPower(const Scalar L, const Scalar m, const Scalar rho)
    {
      const Scalar cond = L / m;
      const Scalar sqrt_L_m = math::sqrt(L * m);
      const Scalar rho_power = math::log(rho / sqrt_L_m) / math::log(cond);
      return rho_power;
    }

  protected:
    Scalar ratio_primal_dual;
    Scalar rho_increment;
  };

  ///
  /// \brief Implementation of LINEAR ADMM update rule.
  template<typename _Scalar>
  struct ADMMLinearUpdateRuleTpl
  {
    typedef _Scalar Scalar;

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

    ADMMLinearUpdateRuleTpl(const Scalar ratio_primal_dual, const Scalar factor)
    : ratio_primal_dual(ratio_primal_dual)
    , increase_factor(factor)
    , decrease_factor(factor)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(factor > Scalar(1), "factor should be greater than one.");
    }

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

  protected:
    Scalar ratio_primal_dual;
    Scalar increase_factor, decrease_factor;
  };

  ///
  /// \brief Implementation of OSQP ADMM update rule.
  template<typename _Scalar>
  struct ADMMOSQPUpdateRuleTpl
  {
    typedef _Scalar Scalar;

    ADMMOSQPUpdateRuleTpl(const Scalar ratio_primal_dual, const Scalar eps_reg)
    : ratio_primal_dual(ratio_primal_dual)
    , eps_reg(eps_reg)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        ratio_primal_dual > Scalar(0), "ratio_primal_dual should be positive.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(eps_reg > Scalar(0), "eps_reg should be positive.");
    }

    void eval(const Scalar primal_feasibility, const Scalar dual_feasibility, Scalar & rho) const
    {
      if (
        primal_feasibility > this->ratio_primal_dual * dual_feasibility //
        || dual_feasibility > this->ratio_primal_dual * primal_feasibility)
      {
        rho *= std::sqrt(primal_feasibility / (dual_feasibility + this->eps_reg));
      }
    }

  protected:
    Scalar ratio_primal_dual;
    Scalar eps_reg;
  };

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

  ///
  /// \brief Settings for the ADMM constraint solver loop.
  template<typename _Scalar>
  struct ADMMSolverSettingsTpl : ConstraintSolverSettingsBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverSettingsBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;

    /// \brief Default constructor
    ADMMSolverSettingsTpl(
      std::size_t max_iterations = 1000,
      Scalar tol_feasibility = Scalar(1e-6),
      Scalar tol_rel_feasibility = Scalar(1e-6),
      Scalar tol_complementarity = Scalar(1e-6),
      Scalar tol_rel_complementarity = Scalar(1e-6),
      bool solve_ncp = true,
      bool measure_timings = false,
      bool stat_record = false,
      std::optional<RefConstVectorXs> preconditioner = std::nullopt,
      std::optional<RefConstVectorXs> primal_guess = std::nullopt,
      std::optional<RefConstVectorXs> dual_guess = std::nullopt,
      std::optional<Scalar> rho_init = std::nullopt,
      ADMMUpdateRule admm_update_rule = ADMMUpdateRule::OSQP,
      ADMMProximalRule admm_proximal_rule = ADMMProximalRule::MANUAL,
      Scalar mu_prox = 1e-6,
      Scalar tau_prox = Scalar(1),
      Scalar tau = Scalar(0.5),
      Scalar ratio_primal_dual = Scalar(10),
      Scalar dual_momentum = Scalar(0),
      Scalar rho_update_ratio = Scalar(0),
      std::size_t rho_min_update_frequency = 1,
      Scalar rho_momentum = Scalar(0),
      Scalar spectral_rho_power_init = Scalar(0.2),
      Scalar spectral_rho_power_factor = Scalar(0.05),
      Scalar linear_update_rule_factor = Scalar(2),
      std::size_t lanczos_size = std::numeric_limits<int>::max(),
      std::size_t max_delassus_decomposition_updates = std::numeric_limits<int>::max(),
      std::size_t anderson_capacity = 0)
    : Base(
        max_iterations,
        tol_feasibility,
        tol_rel_feasibility,
        tol_complementarity,
        tol_rel_complementarity,
        solve_ncp,
        measure_timings,
        stat_record)
    , preconditioner(preconditioner)
    , primal_guess(primal_guess)
    , dual_guess(dual_guess)
    , rho_init(rho_init)
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

    /// \brief Optional preconditioner. Not used yet.
    std::optional<RefConstVectorXs> preconditioner;

    /// \brief Optional guess for the primal variable (impulses).
    std::optional<RefConstVectorXs> primal_guess;

    /// \brief Optional guess for the dual variable (velocities).
    std::optional<RefConstVectorXs> dual_guess;

    /// \brief Initial value of rho parameter.
    /// If set to boost::none, the initial rho will be computed by estimating
    /// the largest eigenvalue of the Delassus.
    std::optional<Scalar> rho_init;

    // ----------------------
    // ADMM specific settings

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

  typedef ADMMSolverSettingsTpl<context::Scalar> ADMMSolverSettings;

  ///
  /// \brief Workspace for the ADMM constraint solver.
  template<typename _Scalar>
  struct ADMMSolverWorkspaceTpl
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
    typedef LanczosDecompositionTpl<MatrixXs> LanczosDecomposition;
    typedef AndersonAccelerationTpl<Scalar> AndersonAcceleration;

    // TODO: use eigenstorage

    static constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();

    /// \brief Default constructor
    ADMMSolverWorkspaceTpl()
    : problem_size(0)
    , lanczos_size(2)
    , anderson_capacity(0)
    , rho(nan)
    , spectral_rho_power(nan)
    , mu_prox(nan)
    , lanczos_decomposition(Eigen::Index(2), Eigen::Index(2))
    , anderson_history(problem_size, anderson_capacity)
    {
    }

    /// \brief Constructor given problem_size, lanczos_size and anderson_capacity.
    ADMMSolverWorkspaceTpl(
      std::size_t problem_size, std::size_t lanczos_size = 2, std::size_t anderson_capacity = 0)
    : problem_size(problem_size)
    , lanczos_size(lanczos_size)
    , anderson_capacity(anderson_capacity)
    , rho(nan)
    , spectral_rho_power(nan)
    , mu_prox(nan)
    , lanczos_decomposition(
        static_cast<Eigen::Index>(math::max(std::size_t(2), problem_size)), //
        static_cast<Eigen::Index>(math::max(std::size_t(2), math::min(problem_size, lanczos_size))))
    , anderson_history(problem_size, anderson_capacity)
    {
      resize(problem_size, lanczos_size, anderson_capacity);
    }

    /// \brief Resize workspace vectors and operators to problem sizes.
    void
    resize(std::size_t problem_size_, std::size_t lanczos_size_, std::size_t anderson_capacity_)
    {
      problem_size = problem_size_;
      lanczos_size = math::max(std::size_t(2), math::min(problem_size, lanczos_size_));
      anderson_capacity = anderson_capacity_;

      const Eigen::Index np = static_cast<Eigen::Index>(problem_size);
      x.setZero(np);
      x_anderson.setZero(np);
      y.setZero(np);
      x_previous.setZero(np);
      y_previous.setZero(np);
      z.setZero(np);
      z_anderson.setZero(np);
      z_previous.setZero(np);
      desaxce.setZero(np);
      rhs.resize(np);
      tmp.resize(np);
      primal_feasibility_vector.setZero(np);
      anderson_primal_feasibility_vector.setZero(np);
      dual_feasibility_vector.setZero(np);

      rho = nan;
      spectral_rho_power = nan;
      mu_prox = nan;

      // resize lanczos
      const std::size_t lanczos_problem_size = math::max(std::size_t(2), problem_size);
      if (
        lanczos_decomposition.size() != static_cast<Eigen::Index>(lanczos_problem_size)
        || lanczos_decomposition.decompositionSize() != static_cast<Eigen::Index>(lanczos_size))
      {
        lanczos_decomposition = LanczosDecomposition(
          static_cast<Eigen::Index>(lanczos_problem_size), static_cast<Eigen::Index>(lanczos_size));
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

    /// \brief Primal variable at current iteration.
    VectorXs x;

    /// \brief Primal variable at previous iteration.
    VectorXs x_previous;

    /// \brief Anderson accelerated primal variable at current iteration.
    VectorXs x_anderson;

    /// \brief Projected primal variable at current iteration.
    VectorXs y;

    /// \brief Projected primal variable at previous iteration.
    VectorXs y_previous;

    /// \brief Dual variable at current iteration.
    VectorXs z;

    /// \brief Dual variable at previous iteration.
    VectorXs z_previous;

    /// \brief Anderson accelerated dual variable at current iteration.
    VectorXs z_anderson;

    /// \brief Desaxce correction (always 0 when solving CCP) at current iteration.
    VectorXs desaxce;

    /// \brief Temporary variable for solving linear system.
    VectorXs rhs;

    /// \brief Temporary variable for holding various vectors.
    VectorXs tmp;

    /// \brief Primal feasibility vector at current iteration.
    VectorXs primal_feasibility_vector;

    /// \brief Anderson related primal feasibility vector at current iteration.
    VectorXs anderson_primal_feasibility_vector;

    /// \brief Dual feasibility vector at current iteration.
    VectorXs dual_feasibility_vector;

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
  }; // struct ADMMSolverWorkspaceTpl

  typedef ADMMSolverWorkspaceTpl<context::Scalar> ADMMSolverWorkspace;

  ///
  /// \brief Struct describing the solution of the ADMM constraint solver
  /// after calling the `solve` method.
  template<typename _Scalar>
  struct ADMMSolverSolutionTpl : ConstraintSolverSolutionBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ConstraintSolverSolutionBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    using Base::nan;

    /// \brief Default constructor.
    ADMMSolverSolutionTpl()
    : Base()
    , problem_size(0)
    , delassus_decomposition_update_count(0)
    , rho(nan)
    , spectral_rho_power(nan)
    , mu_prox(nan)
    {
    }

    /// \brief Reset the solution.
    void reset()
    {
      Base::reset();
      problem_size = 0;
      delassus_decomposition_update_count = 0;
      rho = nan;
      spectral_rho_power = nan;
      mu_prox = nan;
    }

    /// \brief Resize the primal/dual/desaxce vectors of the solution.
    void resize(std::size_t problem_size_)
    {
      problem_size = problem_size_;

      Eigen::Index np = static_cast<Eigen::Index>(problem_size);
      x.setZero(np);
      y.setZero(np);
      z.setZero(np);
      desaxce.setZero(np);
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

    /// \brief Value of ADMM rho term.
    Scalar rho;

    /// \brief Value of ADMM spectral rule rho power at the solution.
    /// This is relevant only if `SPECTRAL` was selected as an update rule.
    Scalar spectral_rho_power;

    /// \brief Value of ADMM proximal term.
    Scalar mu_prox;

    /// \brief Non-projected primal solution.
    VectorXs x;

    /// \brief Primal solution projected onto constraints.
    VectorXs y;

    /// \brief Dual solution.
    VectorXs z;

    /// \brief Desaxce term of the solution
    VectorXs desaxce;
  }; // struct ADMMSolverSolutionTpl

  typedef ADMMSolverSolutionTpl<context::Scalar> ADMMSolverSolution;

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

  typedef ADMMSolverStatsTpl<context::Scalar> ADMMSolverStats;

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

    typedef ADMMSolverWorkspaceTpl<Scalar> ADMMSolverWorkspace;
    typedef ADMMSolverSettingsTpl<Scalar> ADMMSolverSettings;
    typedef ADMMSolverSolutionTpl<Scalar> ADMMSolverSolution;
    typedef ADMMSolverStatsTpl<Scalar> ADMMSolverStats;

    typedef ADMMSpectralUpdateRuleTpl<Scalar> ADMMSpectralUpdateRule;
    typedef ADMMOSQPUpdateRuleTpl<Scalar> ADMMOSQPUpdateRule;
    typedef ADMMLinearUpdateRuleTpl<Scalar> ADMMLinearUpdateRule;
    typedef ADMMUpdateRuleContainerTpl<Scalar> ADMMUpdateRuleContainer;

    /// \brief Default constructor.
    explicit ADMMConstraintSolverTpl(std::size_t problem_size)
    : Base()
    , solution()
    , workspace(problem_size)
    , stats()
    {
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// TODO: put correct params
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] preconditioner Precondtionner of the problem.
    /// \param[in] primal_guess Optional initial guess of the primal solution (constrained forces).
    /// \param[in] dual_guess Optinal Initial guess of the dual solution (constrained velocities).
    /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false)
    /// \param[in] admm_update_rule update rule for ADMM (constant, linear or spectral)
    /// \param[in] rho0 Initial value of the rho parameter.
    /// \param[in] stat_record record solver metrics
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
      const ADMMSolverSettings & settings);

#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL

    ADMMSolverSolution solution;
    ADMMSolverWorkspace workspace;
    ADMMSolverStats stats;

  protected:
    /// \brief Compute largest eigen value of delassus.
    template<typename DelassusDerived>
    static Scalar computeDelassusLargestEigenvalue(
      const DelassusOperatorBase<DelassusDerived> & delassus, ADMMSolverWorkspace & workspace);

    /// \brief Retrieve primal and/or dual guesses from settings.
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
      ADMMSolverWorkspace & workspace,
      const ADMMSolverSettings & settings);
  }; // struct ADMMConstraintSolverTpl

} // namespace pinocchio

#include "pinocchio/algorithm/admm-solver.hxx"

#endif // ifndef __pinocchio_algorithm_admm_solver_hpp__
