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

#include "pinocchio/math/lanczos-decomposition.hpp"

#include "pinocchio/algorithm/diagonal-preconditioner.hpp"

#include <boost/optional.hpp>

namespace pinocchio
{
  template<typename Scalar>
  struct ADMMContactSolverTpl;
  typedef ADMMContactSolverTpl<context::Scalar> ADMMContactSolver;

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
  };

  template<typename _Scalar>
  struct AndersonAccelerationTpl
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    AndersonAccelerationTpl(int problem_size, std::size_t capacity)
    : details(problem_size, capacity)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(capacity >= 0, "capacity needs to be positive");
      this->reserve(problem_size, capacity);
    }

    /// \brief Reserve the capacity of the Anderson acceleration
    void reserve(int new_problem_size, std::size_t new_capacity)
    {
      this->details.problem_size = new_problem_size;
      this->details.capacity = new_capacity;
      this->details.xs.resize(this->capacity(), VectorXs::Zero(this->problem_size()));
      this->details.zs.resize(this->capacity(), VectorXs::Zero(this->problem_size()));
      this->details.zdiffs.resize(this->capacity(), VectorXs::Zero(this->problem_size()));
      this->details.weights.resize(math::max(0, int(this->capacity() - 1)));
      this->details.M.resize(this->problem_size(), math::max(0, int(this->capacity() - 1)));
      this->clear();
    }

    /// \brief Clear the Anderson acceleration.
    void clear()
    {
      this->details.size = 0;
      this->details.idx = 0;
    }

    /// \brief Get the current Anderson acceleration size.
    std::size_t size() const
    {
      return this->details.size;
    }

    /// \brief Get the capacity of this Anderson acceleration.
    std::size_t capacity() const
    {
      return this->details.capacity;
    }

    /// \brief Get the problem size which this Anderson acceleration fits.
    int problem_size() const
    {
      return this->details.problem_size;
    }

    /// \brief Getter for the Anderson weights
    Eigen::VectorBlock<VectorXs> weights()
    {
      return this->details.weights.head(this->size() - 1);
    }

    /// \brief Const getter for the Anderson weights
    Eigen::VectorBlock<const VectorXs> weights() const
    {
      return this->details.weights.head(this->size() - 1);
    }

    /// \brief Push back default iterates into the Anderson acceleration.
    template<typename VectorLikeX, typename VectorLikeZ, typename VectorLikeZDiff>
    void push_back(
      const Eigen::MatrixBase<VectorLikeX> & x,
      const Eigen::MatrixBase<VectorLikeZ> & z,
      const Eigen::MatrixBase<VectorLikeZDiff> & zdiff)
    {
      if (this->size() > 0)
      {
        // cycle through the std::vector to maintain correct history
        this->details.idx = (this->details.idx + 1) % this->capacity();
      }
      else
      {
        assert(this->details.idx == 0);
      }

      const std::size_t idx = this->details.idx;
      this->details.xs[idx] = x;
      this->details.zs[idx] = z;
      this->details.zdiffs[idx] = zdiff;

      // update Anderson current history size
      this->details.size = math::min(this->size() + 1, this->capacity());
    }

    /// \brief Fit the Anderson acceleration and store results in weights.
    void fit()
    {
      if (this->size() < 2 || this->size() < this->capacity())
        return;

      const std::size_t idx = this->details.idx;
      const auto & zdiffs = this->details.zdiffs;
      auto M = this->details.M.leftCols(this->size() - 1);

      for (std::size_t i = 0; i < this->size() - 1; ++i)
      {
        std::size_t i1 = (idx - i) % this->capacity();
        std::size_t i2 = (idx - i - 1) % this->capacity();
        M.col(int(i)) = zdiffs[i1] - zdiffs[i2];
      }

      // fit the Anderson weights
      this->weights() = M.colPivHouseholderQr().solve(zdiffs[idx]);
    }

    template<typename VectorLikeX, typename VectorLikeZ>
    void getAcceleratedIterates(
      const Eigen::MatrixBase<VectorLikeX> & x_, const Eigen::MatrixBase<VectorLikeZ> & z_) const
    {
      VectorLikeX & x = x_.const_cast_derived();
      VectorLikeZ & z = z_.const_cast_derived();

      const std::size_t idx = this->details.idx;
      const auto & xs = this->details.xs;
      const auto & zs = this->details.zs;

      x = xs[idx];
      z = zs[idx];

      if (this->size() < 2 || this->size() < this->capacity())
        return;

      for (std::size_t i = 0; i < this->size() - 1; ++i)
      {
        const std::size_t i1 = (idx - i) % this->capacity();
        const std::size_t i2 = (idx - i - 1) % this->capacity();
        x -= this->weights().coeff(int(i)) * (xs[i1] - xs[i2]);
        z -= this->weights().coeff(int(i)) * (zs[i1] - zs[i2]);
      }
    }

    struct AndersonAccelerationData
    {
      std::vector<VectorXs> xs;     // history of x (first primal variable)
      std::vector<VectorXs> zs;     // history of z (dual variable)
      std::vector<VectorXs> zdiffs; // history of dual residuals
      VectorXs weights;             // weights of Anderson acceleration, computed by `fit`
      MatrixXs M;                   // matrix used to fit Anderson acceleration weights
      int problem_size;             // size of each history vector
      std::size_t capacity;         // capacity of the history
      std::size_t size;             // size of the history
      std::size_t idx;              // index of the most recent element in history

      AndersonAccelerationData(int problem_size, std::size_t capacity)
      : problem_size(problem_size)
      , capacity(capacity)
      , size(0)
      , idx(0)
      {
      }
    };

    /// \brief Internal details for Anderson acceleration.
    /// An experienced external user should only read this data e.g. for debug purposes.
    AndersonAccelerationData details;
  };

  template<typename _Scalar>
  struct PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
    ADMMContactSolverTpl : ContactSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ContactSolverBaseTpl<_Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Ref<VectorXs> RefVectorXs;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;
    typedef const Eigen::Ref<const VectorXs> ConstRefVectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
    typedef LanczosDecompositionTpl<MatrixXs> LanczosDecomposition;
    typedef DiagonalPreconditionerTpl<VectorXs> DiagonalPreconditioner;
    typedef AndersonAccelerationTpl<Scalar> AndersonAcceleration;

    using Base::problem_size;

    //    struct SolverParameters
    //    {
    //      explicit SolverParameters(const int problem_dim)
    //      : rho_power(Scalar(0.2))
    //      , ratio_primal_dual(Scalar(10))
    //      , mu_prox
    //      {
    //
    //      }
    //
    //      /// \brief Rho solver ADMM
    //      boost::optional<Scalar> rho;
    //      /// \brief Power value associated to rho. This quantity will be automatically updated.
    //      Scalar rho_power;
    //      /// \brief Ratio primal/dual
    //      Scalar ratio_primal_dual;
    //      /// \brief Proximal value
    //      Scalar mu_prox;
    //
    //      /// \brief Largest eigenvalue
    //      boost::optional<Scalar> L_value;
    //      /// \brief Largest eigenvector
    //      boost::optional<VectorXs> L_vector;
    //    };
    //
    struct ADMMSolverStats : Base::SolverStats
    {
      ADMMSolverStats()
      : Base::SolverStats()
      , delassus_decomposition_update_count(0)
      {
      }

      explicit ADMMSolverStats(const int max_it)
      : Base::SolverStats(max_it)
      , delassus_decomposition_update_count(0)
      {
        reserve(max_it);
      }

      void reserve(const int max_it)
      {
        rho.reserve(size_t(max_it));
        mu_prox.reserve(size_t(max_it));
        anderson_size.reserve(size_t(max_it));
        linear_system_residual.reserve(size_t(max_it));
        linear_system_consistency.reserve(size_t(max_it));
      }

      void reset()
      {
        Base::SolverStats::reset();
        rho.clear();
        mu_prox.clear();
        anderson_size.clear();
        linear_system_residual.clear();
        linear_system_consistency.clear();
        delassus_decomposition_update_count = 0;
      }

      ///  \brief Number of Delassus decomposition updates.
      int delassus_decomposition_update_count;

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
    };

    //
    //    struct SolverResults
    //    {
    //      explicit SolverResults(const int problem_dim, const int max_it)
    //      : L_vector(problem_dim)
    //
    //      /// \brief Largest eigenvalue
    //      Scalar L_value;
    //      /// \brief Largest eigenvector
    //      VectorXs L_vector;
    //
    //      SolverStats stats;
    //    };

    explicit ADMMContactSolverTpl(
      int problem_dim,
      Scalar tau_prox = Scalar(1),
      Scalar tau = Scalar(0.5),
      Scalar rho_power = Scalar(0.2),
      Scalar rho_power_factor = Scalar(0.05),
      Scalar linear_update_rule_factor = Scalar(2),
      Scalar ratio_primal_dual = Scalar(10),
      int lanczos_size = int(3),
      int max_delassus_decomposition_updates = std::numeric_limits<int>::max(),
      Scalar dual_momentum = Scalar(0),
      Scalar rho_momentum = Scalar(0),
      Scalar rho_update_ratio = Scalar(0),
      int rho_min_update_frequency = int(1),
      std::size_t anderson_capacity = std::size_t(0))
    : Base(problem_dim)
    , is_initialized(false)
    , tau_prox(tau_prox)
    , mu_prox(Scalar(1e-6))
    , tau(tau)
    , rho(Scalar(10))
    , rho_power(rho_power)
    , rho_power_factor(rho_power_factor)
    , linear_update_rule_factor(linear_update_rule_factor)
    , ratio_primal_dual(ratio_primal_dual)
    , lanczos_decomposition(
        static_cast<Eigen::DenseIndex>(math::max(2, problem_dim)),
        static_cast<Eigen::DenseIndex>(math::max(2, math::min(lanczos_size, problem_dim))))
    , x_(VectorXs::Zero(problem_dim))
    , x_anderson_(VectorXs::Zero(problem_dim))
    , y_(VectorXs::Zero(problem_dim))
    , x_previous_(VectorXs::Zero(problem_dim))
    , y_previous_(VectorXs::Zero(problem_dim))
    , z_(VectorXs::Zero(problem_dim))
    , z_anderson_(VectorXs::Zero(problem_dim))
    , z_previous_(VectorXs::Zero(problem_dim))
    , s_(VectorXs::Zero(problem_dim))
    , rhs(problem_dim)
    , tmp(problem_dim)
    , primal_feasibility_vector(VectorXs::Zero(problem_dim))
    , anderson_primal_feasibility_vector(VectorXs::Zero(problem_dim))
    , dual_feasibility_vector(VectorXs::Zero(problem_dim))
    , delassus_decomposition_update_count(0)
    , max_delassus_decomposition_updates(max_delassus_decomposition_updates)
    , dual_momentum(dual_momentum)
    , rho_momentum(rho_momentum)
    , rho_update_ratio(rho_update_ratio)
    , rho_min_update_frequency(rho_min_update_frequency)
    , anderson_history(problem_size, anderson_capacity)
    , stats()
    {
    }

    /// \brief Reset the solver.
    void reset()
    {
      this->it = 0;
      this->delassus_decomposition_update_count = 0;
      this->stats.reset();
      this->is_initialized = false;
    }

    /// \brief Get the ADMM penalty value.
    Scalar getRho() const
    {
      return this->rho;
    }

    /// \brief Set the power associated to the problem conditionning.
    void setRhoPower(const Scalar rho_power)
    {
      this->rho_power = rho_power;
    }
    /// \brief Get the power associated to the problem conditionning.
    Scalar getRhoPower() const
    {
      return rho_power;
    }

    /// \brief Set the power factor associated to the problem conditionning.
    /// Only related to ADMMUpdateRule::SPECTRAL.
    void setRhoPowerFactor(const Scalar rho_power_factor)
    {
      this->rho_power_factor = rho_power_factor;
    }
    /// \brief Get the value of the increase/decrease factor associated to the problem
    /// conditionning.
    /// Only related to ADMMUpdateRule::SPECTRAL.
    Scalar getRhoPowerFactor() const
    {
      return rho_power_factor;
    }

    /// \brief Set rho momentum (rho = pow(rho, momentum) * pow(new_rho, 1 - momentum)).
    /// Value of 0 is no momentum.
    void setRhoMomentum(const Scalar rho_momentum)
    {
      this->rho_momentum = rho_momentum;
    }
    /// \brief Get rho momentum (rho = pow(rho, momentum) * pow(new_rho, 1 - momentum)).
    /// Value of 0 is no momentum.
    Scalar getRhoMomentum() const
    {
      return this->rho_momentum;
    }

    /// \brief Set the update factor of the Linear update rule
    /// Only related to ADMMUpdateRule::LINEAR.
    void setLinearUpdateRuleFactor(const Scalar linear_update_rule_factor)
    {
      this->linear_update_rule_factor = linear_update_rule_factor;
    }
    /// \brief Get the value of the increase/decrease factor of the Linear update rule
    /// Only related to ADMMUpdateRule::LINEAR.
    Scalar getLinearUpdateRuleFactor() const
    {
      return linear_update_rule_factor;
    }

    /// \brief Set the tau linear scaling factor.
    void setTau(const Scalar tau)
    {
      this->tau = tau;
    }
    /// \brief Get the tau linear scaling factor.
    Scalar getTau() const
    {
      return tau;
    }

    /// \brief Set the tau linear proximal scaling factor.
    void setProximalTau(const Scalar tau_prox)
    {
      this->tau_prox = tau_prox;
    }
    /// \brief Get the tau linear proximal scaling factor.
    Scalar getProximalTau() const
    {
      return this->tau_prox;
    }

    /// \brief Get the proximal value.
    Scalar getProximalValue() const
    {
      return this->mu_prox;
    }

    /// \brief Set the maximum number of decompositions of the Delassus.
    void setMaxDelassusDecompositionUpdates(const int max_delassus_decomposition_updates)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        max_delassus_decomposition_updates >= int(0),
        "max_delassus_decomposition_updates should be positive");
      this->max_delassus_decomposition_updates = max_delassus_decomposition_updates;
    }
    /// \brief Get the maximum number of decompositions of the Delassus.
    int getMaxDelassusDecompositionUpdates() const
    {
      return this->max_delassus_decomposition_updates;
    }

    /// \brief Set dual momentum.
    void setDualMomentum(const Scalar dual_momentum)
    {
      this->dual_momentum = dual_momentum;
    }
    /// \brief Get dual momentum.
    Scalar getDualMomentum() const
    {
      return this->dual_momentum;
    }

    /// \brief Set the rho update ratio.
    /// The rho is only updated if the ratio between the current rho and the new one is bigger/lower
    /// than a threshold ratio.
    void setRhoUpdateRatio(const Scalar rho_update_ratio)
    {
      this->rho_update_ratio = rho_update_ratio;
    }
    /// \brief Get the rho update ratio.
    Scalar getRhoUpdateRatio() const
    {
      return this->rho_update_ratio;
    }

    /// \brief Set the minimum rho update frequency.
    /// \copydoc rho_min_update_frequency
    void setRhoMinUpdateFrequency(const Scalar rho_min_update_frequency)
    {
      this->rho_min_update_frequency = rho_min_update_frequency;
    }
    /// \brief Get the minimum rho update frequency.
    /// \copydoc rho_min_update_frequency
    Scalar getRhoMinUpdateFrequency() const
    {
      return this->rho_min_update_frequency;
    }

    /// \brief Set the primal/dual ratio.
    void setRatioPrimalDual(const Scalar ratio_primal_dual)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        ratio_primal_dual > 0., "The ratio primal/dual should be positive strictly");
      this->ratio_primal_dual = ratio_primal_dual;
    }
    /// \brief Get the primal/dual ratio.
    Scalar getRatioPrimalDual() const
    {
      return ratio_primal_dual;
    }

    ///  \returns the number of updates of the Delassus decomposition.
    int getDelassusDecompositionUpdateCount() const
    {
      return this->delassus_decomposition_update_count;
    }

    /// \brief Sets the size of triangular matrix of Lanczos decomposition.
    /// The higher the size, the more accurate the estimation of min/max eigenvalues will be.
    /// Note: the maximum size is the size of the problem
    void setLanczosSize(int decomposition_size)
    {
      // TODO(lmontaut): should we throw if size > problem_size or instead take the min as done
      // below?
      int new_lanczos_size = math::max(2, this->problem_size);
      int new_lanczos_decomposition_size =
        math::max(2, math::min(decomposition_size, this->problem_size));
      if (
        new_lanczos_size != this->lanczos_decomposition.size()
        || new_lanczos_decomposition_size != this->lanczos_decomposition.decompositionSize())
      {
        this->lanczos_decomposition = LanczosDecomposition(
          static_cast<Eigen::DenseIndex>(new_lanczos_size),
          static_cast<Eigen::DenseIndex>(new_lanczos_decomposition_size));
      }
    }

    /// \returns the Lanczos algorithm used for eigenvalues estimation.
    const LanczosDecomposition & getLanczosDecomposition() const
    {
      return lanczos_decomposition;
    }

    /// \brief Get the capacity of the Anderson acceleration.
    /// \copydoc anderson_history
    std::size_t getAndersonAccelerationCapacity() const
    {
      return this->anderson_history.capacity();
    }
    ///
    /// \brief Set the capacity of the Anderson acceleration.
    /// \copydoc anderson_history
    void setAndersonAccelerationCapacity(const std::size_t anderson_history_capacity)
    {
      return this->anderson_history.reserve(this->problem_size, anderson_history_capacity);
    }

    ADMMSolverStats & getStats()
    {
      return stats;
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
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
      template<typename T> class Holder,
      typename ConstraintModel,
      typename ConstraintAllocator>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<Holder<const ConstraintModel>, ConstraintAllocator> & constraint_models,
      const Scalar dt,
      const boost::optional<RefConstVectorXs> preconditioner = boost::none,
      const boost::optional<RefConstVectorXs> primal_guess = boost::none,
      const boost::optional<RefConstVectorXs> dual_guess = boost::none,
      const bool solve_ncp = true,
      const ADMMUpdateRule admm_update_rule = ADMMUpdateRule::SPECTRAL,
      const boost::optional<Scalar> rho0 = boost::none,
      const ADMMProximalRule admm_proximal_policy = ADMMProximalRule::MANUAL,
      const boost::optional<Scalar> mu_prox0 = boost::none,
      const bool stat_record = false);

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
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
      typename ConstraintAllocator>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintAllocator> & constraint_models,
      const Scalar dt,
      const boost::optional<RefConstVectorXs> preconditioner = boost::none,
      const boost::optional<RefConstVectorXs> primal_guess = boost::none,
      const boost::optional<RefConstVectorXs> dual_guess = boost::none,
      const bool solve_ncp = true,
      const ADMMUpdateRule admm_update_rule = ADMMUpdateRule::SPECTRAL,
      const boost::optional<Scalar> rho0 = boost::none,
      const ADMMProximalRule admm_proximal_policy = ADMMProximalRule::MANUAL,
      const boost::optional<Scalar> mu_prox0 = boost::none,
      const bool stat_record = false)
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      return solve(
        delassus, g, wrapped_constraint_models, dt, preconditioner, primal_guess, dual_guess,
        solve_ncp, admm_update_rule, rho0, admm_proximal_policy, mu_prox0, stat_record);
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] primal_guess Optional initial guess of the primal solution (constrained forces).
    /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false)
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      template<typename T> class Holder,
      typename ConstraintModel,
      typename ConstraintAllocator,
      typename VectorLikeOut>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<Holder<const ConstraintModel>, ConstraintAllocator> & constraint_models,
      const Scalar dt,
      const Eigen::DenseBase<VectorLikeOut> & primal_guess,
      const bool solve_ncp = true)
    {
      return solve(
        delassus.derived(), g.derived(), constraint_models, dt, boost::none, primal_guess.derived(),
        boost::none, solve_ncp);
    }

    ///
    /// \brief Solve the constraint problem composed of problem data (G,g,constraint_models) and
    /// starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the constraint problem.
    /// \param[in] g Free constraint acceleration or velicity associted with the constraint problem.
    /// \param[in] constraint_models Vector of constraints.
    /// \param[in] primal_guess Optional initial guess of the primal solution (constrained forces).
    /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false)
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLike,
      typename ConstraintModel,
      typename ConstraintAllocator,
      typename VectorLikeOut>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLike> & g,
      const std::vector<ConstraintModel, ConstraintAllocator> & constraint_models,
      const Scalar dt,
      const Eigen::DenseBase<VectorLikeOut> & primal_guess,
      const bool solve_ncp = true)
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        constraint_models.cbegin(), constraint_models.cend());

      return solve(delassus, g, wrapped_constraint_models, dt, primal_guess, solve_ncp);
    }

    /// \returns the primal solution of the problem
    const VectorXs & getPrimalSolution() const
    {
      return y_;
    }
    /// \returns the dual solution of the problem
    const VectorXs & getDualSolution() const
    {
      return z_;
    }
    /// \returns the complementarity shift
    const VectorXs & getComplementarityShift() const
    {
      return s_;
    }

    bool isInitialized() const
    {
      return this->is_initialized;
    }

    template<typename DelassusDerived>
    Scalar computeDelassusLargestEigenvalue(const DelassusOperatorBase<DelassusDerived> & delassus);

  protected:
    /// \brief Default value of ADMM proximal term.
    static constexpr Scalar mu_prox_default = 1e-6;

    bool is_initialized;

    /// \brief Linear scaling of the ADMM proximal term.
    Scalar tau_prox;
    /// \brief Value of the ADMM proximal term.
    Scalar mu_prox;

    /// \brief Linear scaling of the ADMM penalty term
    Scalar tau;
    /// \brief Penalty term associated to the ADMM.
    Scalar rho;

    // Set of parameters associated with the Spectral update rule
    /// \brief Power value associated to rho. This quantity will be automatically updated.
    Scalar rho_power;
    /// \brief Update factor for the primal/dual update of rho.
    Scalar rho_power_factor;

    // Set of parameters associated with the Linear update rule
    /// \brief value of the increase/decrease factor
    Scalar linear_update_rule_factor;

    ///  \brief Ratio primal/dual
    Scalar ratio_primal_dual;

    /// \brief Lanczos decomposition algorithm.
    LanczosDecomposition lanczos_decomposition;

    /// \brief Primal variables (corresponds to the constraint impulses)
    VectorXs x_, x_anderson_, y_, x_previous_, y_previous_;
    /// \brief Dual variable of the ADMM (corresponds to the contact velocity or acceleration).
    VectorXs z_, z_anderson_, z_previous_;
    /// \brief De Saxé shift
    VectorXs s_;

    VectorXs rhs, tmp, primal_feasibility_vector, anderson_primal_feasibility_vector,
      dual_feasibility_vector;

    int delassus_decomposition_update_count;
    int max_delassus_decomposition_updates;

    /// \brief Momentum on the dual
    Scalar dual_momentum;

    /// \brief Momentum of rho (rho = pow(rho, momentum) * pow(new_rho, 1 - momentum)).
    /// Value of 0 is no momentum.
    Scalar rho_momentum;

    /// \brief The rho is only updated if the ratio between the current rho and the new one is
    /// bigger/lower than a threshold ratio.
    Scalar rho_update_ratio;

    /// \brief Rho min update frequency: the solver has to wait at least rho_min_update_frequency
    /// until it can trigger a new rho update.
    Scalar rho_min_update_frequency;

    /// \brief Anderson acceleration history.
    /// An Anderson acceleration of capacity <= 1 is inactive (it is the standard ADMM algorithm).
    /// The Anderson acceleration only triggers if the capacity (and the current Anderson size) is
    /// >= 2.
    AndersonAcceleration anderson_history;

    /// \brief Stats recorded by the solver if `solve` is called with `stat_record = true`.
    ADMMSolverStats stats;

#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL
  }; // struct ADMMContactSolverTpl

} // namespace pinocchio

#include "pinocchio/algorithm/admm-solver.hxx"

#endif // ifndef __pinocchio_algorithm_admm_solver_hpp__
