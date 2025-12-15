//
// Copyright (c) 2022-2025 INRIA
//

#ifndef __pinocchio_algorithm_solvers_admm_solver_hxx__
#define __pinocchio_algorithm_solvers_admm_solver_hxx__

#include <limits>

#include "pinocchio/algorithm/solvers/constraint-solver-utils.hpp"
#include "pinocchio/algorithm/constraints/sets/sets.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"
#include "pinocchio/algorithm/delassus-operator-preconditioned.hpp"
#include "pinocchio/utils/reference.hpp"

#include "pinocchio/tracy.hpp"

#include <Eigen/Eigenvalues>

namespace pinocchio
{
  //
  // template<typename DriftVectorLike, typename Scalar>
  // struct ZeroInitialGuessMaxConstraintViolationVisitor
  // : visitors::ConstraintUnaryVisitorBase<
  //     ZeroInitialGuessMaxConstraintViolationVisitor<DriftVectorLike, Scalar>>
  // {
  //   using ArgsType = boost::fusion::vector<const DriftVectorLike &, Scalar &>;
  //   using Base = visitors::ConstraintUnaryVisitorBase<
  //     ZeroInitialGuessMaxConstraintViolationVisitor<DriftVectorLike, Scalar>>;
  //
  //   template<typename ConstraintModel>
  //   static void algo(
  //     const ConstraintModelBase<ConstraintModel> & cmodel,
  //     const DriftVectorLike & drift,
  //     Scalar & max_violation)
  //   {
  //     return algo_impl(cmodel.set(), drift, max_violation);
  //   }
  //
  //   template<typename VectorLike>
  //   static void algo_impl(
  //     const CoulombFrictionConeTpl<Scalar> & set,
  //     const Eigen::MatrixBase<VectorLike> & drift,
  //     Scalar & max_violation)
  //   {
  //     PINOCCHIO_UNUSED_VARIABLE(set);
  //     const Scalar violation = -drift.coeff(2);
  //     if (violation > max_violation)
  //     {
  //       max_violation = violation;
  //     }
  //   }
  //
  //   template<typename VectorLike>
  //   static void algo_impl(
  //     const FullSpaceConeTpl<Scalar> & set,
  //     const Eigen::MatrixBase<VectorLike> & drift,
  //     Scalar & max_violation)
  //   {
  //     PINOCCHIO_UNUSED_VARIABLE(set);
  //     const Scalar violation = drift.template lpNorm<Eigen::Infinity>();
  //     if (violation > max_violation)
  //     {
  //       max_violation = violation;
  //     }
  //   }
  //
  //   template<typename VectorLike>
  //   static void algo_impl(
  //     const BoxSetTpl<Scalar> & set,
  //     const Eigen::MatrixBase<VectorLike> & drift,
  //     Scalar & max_violation)
  //   {
  //     PINOCCHIO_UNUSED_VARIABLE(set);
  //     const Scalar violation = drift.template lpNorm<Eigen::Infinity>();
  //     if (violation > max_violation)
  //     {
  //       max_violation = violation;
  //     }
  //   }
  //
  //   template<typename VectorLike>
  //   static void algo_impl(
  //     const NonNegativeOrthantConeTpl<Scalar> & set,
  //     const Eigen::MatrixBase<VectorLike> & drift,
  //     Scalar & max_violation)
  //   {
  //     Scalar violation = 0;
  //     if (set.size() > 0)
  //     {
  //       violation = -math::min(Scalar(0), drift.minCoeff());
  //     }
  //
  //     if (violation > max_violation)
  //     {
  //       max_violation = violation;
  //     }
  //   }
  //
  //   template<typename ConstraintSet, typename VectorLike>
  //   static void algo_impl(
  //     const ConstraintSet & set,
  //     const Eigen::MatrixBase<VectorLike> & drift,
  //     Scalar & max_violation)
  //   {
  //     // do nothing
  //     PINOCCHIO_UNUSED_VARIABLE(set);
  //     PINOCCHIO_UNUSED_VARIABLE(drift);
  //     PINOCCHIO_UNUSED_VARIABLE(max_violation);
  //   }
  //
  //   /// ::run for individual constraints
  //   template<typename ConstraintModel>
  //   static void run(
  //     const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
  //     const DriftVectorLike & drift,
  //     Scalar & max_violation)
  //   {
  //     algo(cmodel.derived(), drift, max_violation);
  //   }
  //
  //   /// ::run for constraints variant
  //   template<int Options, template<typename S, int O> class ConstraintCollectionTpl>
  //   static void run(
  //     const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
  //     const DriftVectorLike & drift,
  //     Scalar & max_violation)
  //   {
  //     ArgsType args(drift, max_violation);
  //     // Note: Base::run will call `algo` of this visitor
  //     Base::run(cmodel.derived(), args);
  //   }
  // }; // struct ZeroInitialGuessMaxConstraintViolationVisitor
  //
  // template<
  //   typename ConstraintModel,
  //   typename ConstraintModelAllocator,
  //   typename ConstraintData,
  //   typename ConstraintDataAllocator,
  //   typename VectorLikeIn>
  // typename ConstraintModel::Scalar computeZeroInitialGuessMaxConstraintViolation(
  //   const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
  //   const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
  //   const Eigen::DenseBase<VectorLikeIn> & drift)
  // {
  //   PINOCCHIO_TRACY_ZONE_SCOPED_N("computeZeroInitialGuessMaxConstraintViolation");
  //   assert(
  //     constraint_models.size() == constraint_datas.size()
  //     && "Both std::vector should be of equal size.");
  //
  //   Eigen::DenseIndex cindex = 0;
  //
  //   using SegmentType = typename VectorLikeIn::ConstSegmentReturnType;
  //   using Scalar = typename ConstraintModel::Scalar;
  //
  //   Scalar max_violation = Scalar(0);
  //   for (size_t k = 0; k < constraint_models.size(); ++k)
  //   {
  //     const auto & cmodel = helper::get_ref(constraint_models[k]);
  //     const auto & cdata = helper::get_ref(constraint_datas[k]);
  //     const auto csize = cmodel.residualSize(cdata);
  //
  //     SegmentType drift_segment = drift.segment(cindex, csize);
  //     typedef ZeroInitialGuessMaxConstraintViolationVisitor<SegmentType, Scalar> Algo;
  //
  //     Algo::run(cmodel, drift_segment, max_violation);
  //
  //     cindex += csize;
  //   }
  //   return max_violation;
  // }

  template<typename Scalar>
  template<typename DelassusDerived>
  Scalar ADMMConstraintSolverTpl<Scalar>::computeDelassusLargestEigenvalue(
    const DelassusOperatorBase<DelassusDerived> & delassus, ADMMSolverWorkspace & workspace)
  {
    const DelassusDerived & G = delassus.derived();
    Scalar L = Scalar(-1);
    if (workspace.problem_size > 1)
    {
      PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMConstraintSolverTpl::solve - lanczos");
      workspace.lanczos_decomposition.compute(G);
      L = ::pinocchio::computeLargestEigenvalue(workspace.lanczos_decomposition.Ts(), 1e-8);
#ifndef NDEBUG
      const bool enforce_symmetry = true;
      MatrixXs delassus = G.matrix(enforce_symmetry);
      Eigen::SelfAdjointEigenSolver<MatrixXs> solver(delassus);
      VectorXs eigvals = solver.eigenvalues();
      Scalar true_L = eigvals.maxCoeff();
      PINOCCHIO_UNUSED_VARIABLE(true_L);
        //          if (true_m > 0)
        //          {
        //            assert(
        //              math::fabs((true_m - m) / math::max(true_m, m)) < 0.01
        //              && "true_m and m are too far apart.");
        //          }
        // assert(
        //   math::fabs((true_L - L) / math::max(true_L, L)) < 0.01
        //   && "true_L and L are too far apart.");
#endif // NDEBUG
    }
    else
    {
      // TODO adapt this with Gbar
      typedef Eigen::Matrix<Scalar, 1, 1> Vector1;
      const Vector1 Gvec = G * Vector1::Constant(1);
      L = Gvec.coeff(0);
    }

    assert(L > 0 && "L must be positive.");
    return L;
  }

  template<typename _Scalar>
  template<
    typename DelassusDerived,
    typename VectorLike,
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  bool ADMMConstraintSolverTpl<_Scalar>::solve(
    DelassusOperatorBase<DelassusDerived> & delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const ADMMSolverSettings & settings)
  {
    // for easier access
    static constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();
    ADMMSolverSolution & sol = solution;
    ADMMSolverWorkspace & wk = workspace;
    DelassusDerived & G = delassus.derived();

    // Configure/reset solution, workspace and stats
    const Eigen::Index np = G.rows();
    const std::size_t problem_size = static_cast<std::size_t>(np);
    assert(G.cols() == np);
    assert(g.size() == np);
    assert(residualSize(constraint_models, constraint_datas) == np);
    const Scalar min_compliance = G.getCompliance().minCoeff();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      min_compliance >= Scalar(0), "compliance should be a positive vector.");
    //
    sol.reset();
    sol.resize(problem_size);
    assert(sol.iterations == 0);
    //
    wk.resize(problem_size, settings.lanczos_size, settings.anderson_capacity);
    assert(wk.problem_size == problem_size);
    assert(wk.x.size() == np);
    //
    stats.reset();
    if (settings.stat_record)
    {
      stats.reserve(settings.max_iterations);
    }
    //
    settings.checkValidity();

#ifdef PINOCCHIO_WITH_HPP_FCL
    if (settings.measure_timings)
    {
      timer.start();
    }
#endif // PINOCCHIO_WITH_HPP_FCL

    // First, we initialize some utils
    // Initialize De Saxé shift to 0
    // For the CCP, there is no shift.
    // For the NCP, the shift will be initialized using z.
    wk.desaxce.setZero();

    // Set initial damping of the delassus to the proximal value and get smallest possible
    // eigenvalue of the problem.
    wk.mu_prox = settings.mu_prox;
    Scalar m = min_compliance + wk.mu_prox;
    wk.rhs = VectorXs::Constant(np, wk.mu_prox);
    G.updateDamping(wk.rhs);
    sol.delassus_decomposition_update_count++;

    // Initialization of the primal/dual variables.
    // If both primal and dual guesses are given, the solver uses both.
    // If one is given but not the other, the solver will compute the missing one using the given
    // one.
    retrievePrimalDualGuess(G, g, constraint_models, constraint_datas, wk, settings);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.x.size(), np);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.y.size(), np);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.z.size(), np);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.desaxce.size(), np);

    // Check NCP/CCP conditions. If they are satisfied, don't run the solver.
    // -- always primaly feasible as y is projected onto the constraints
    sol.primal_feasibility = Scalar(0);
    // -- dual feasibility
    G.applyOnTheRight(wk.y, wk.rhs);
    wk.rhs += g - wk.y.cwiseProduct(G.getDamping());
    if (settings.solve_ncp)
    {
      internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.rhs, wk.desaxce);
      wk.rhs += wk.desaxce;
    }
    internal::computeDualConeProjection(constraint_models, constraint_datas, wk.rhs, wk.tmp);
    wk.tmp -= wk.rhs;
    sol.dual_feasibility = wk.tmp.template lpNorm<Eigen::Infinity>();
    // -- complementarity
    internal::computeConicComplementarity(
      constraint_models, constraint_datas, wk.rhs, wk.y, sol.complementarity);

    bool abs_prec_reached = false;
    bool rel_prec_reached = false;
    if (
      check_expression_if_real<Scalar, false>(sol.complementarity <= settings.tol_complementarity)
      && check_expression_if_real<Scalar, false>(sol.dual_feasibility <= settings.tol_feasibility))
    {
      abs_prec_reached = true;
      wk.z = wk.rhs; // store dual solution
    }

    // init rho power of spectral rule
    wk.spectral_rho_power = settings.spectral_rho_power_init;

    // init rho
    Scalar L = Scalar(-1); // not yet computed
    bool delassus_largest_eigenvalue_computed = false;
    if (settings.rho_init)
    {
      wk.rho = settings.rho_init.value();
    }
    else
    {
      // Compute rho with spectral rule
      L = computeDelassusLargestEigenvalue(G, wk);
      delassus_largest_eigenvalue_computed = true;
      if (std::isnan(L))
      {
        L = Scalar(1);
      }
      wk.rho = ADMMSpectralUpdateRule::computeRho(L, m, wk.spectral_rho_power);
    }
    PINOCCHIO_CHECK_INPUT_ARGUMENT(wk.rho >= 0, "rho should be positive.");
    // clamp the rho
    const Scalar rho_min = 1e-6;
    const Scalar rho_max = 1e6;
    wk.rho = math::max(math::min(wk.rho, rho_max), rho_min);

    // set mu_prox according to prox policy
    switch (settings.admm_proximal_rule)
    {
    case (ADMMProximalRule::MANUAL):
      wk.mu_prox = settings.mu_prox;
      break;
    case (ADMMProximalRule::AUTOMATIC):
      wk.mu_prox = wk.rho;
      break;
    }
    PINOCCHIO_CHECK_INPUT_ARGUMENT(wk.mu_prox >= 0, "mu_prox should be positive.");

    if (!abs_prec_reached)
    {
      // Setup ADMM update rules:
      // Before running ADMM, we compute the largest and smallest eigenvalues of delassus in order
      // to be able to use a spectral update rule for the proximal parameter (rho)
      ADMMUpdateRuleContainer admm_update_rule_container;
      switch (settings.admm_update_rule)
      {
      case (ADMMUpdateRule::SPECTRAL): {
        if (!delassus_largest_eigenvalue_computed)
        {
          // L has not yet been computed, we compute it
          L = computeDelassusLargestEigenvalue(G, wk);
          delassus_largest_eigenvalue_computed = true;
        }
        admm_update_rule_container.spectral_rule = ADMMSpectralUpdateRule(
          settings.ratio_primal_dual, L, m, settings.spectral_rho_power_factor);
        break;
      }
      case (ADMMUpdateRule::OSQP):
        admm_update_rule_container.osqp_rule = ADMMOSQPUpdateRule(settings.ratio_primal_dual, 1e-8);
        break;
      case (ADMMUpdateRule::LINEAR):
        admm_update_rule_container.linear_rule =
          ADMMLinearUpdateRule(settings.ratio_primal_dual, settings.linear_update_rule_factor);
        break;
      case (ADMMUpdateRule::CONSTANT):
        break;
      }

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();

      // Update the decomposition of the Delassus
      Scalar prox_value = settings.tau_prox * wk.mu_prox + settings.tau * wk.rho;
      wk.rhs = VectorXs::Constant(np, prox_value);
      G.updateDamping(wk.rhs);
      Scalar old_prox_value = prox_value;
      sol.delassus_decomposition_update_count++;

      // End of Initialization phase
      wk.x_anderson = wk.x;
      wk.z_anderson = wk.z;
      wk.anderson_history.clear();
      Scalar anderson_primal_feasibility;
      Scalar anderson_previous_primal_feasibility = std::numeric_limits<Scalar>::max();

      Scalar dx_norm = nan;
      Scalar dy_norm = nan;
      Scalar dz_norm = nan;
      const Scalar g_norm_inf = g.template lpNorm<Eigen::Infinity>();
      Scalar x_norm_inf = wk.x.template lpNorm<Eigen::Infinity>();
      Scalar y_norm_inf = wk.y.template lpNorm<Eigen::Infinity>();
      Scalar z_norm_inf = wk.z.template lpNorm<Eigen::Infinity>();
      Scalar x_previous_norm_inf = x_norm_inf;
      Scalar y_previous_norm_inf = y_norm_inf;
      Scalar z_previous_norm_inf = z_norm_inf;

      sol.iterations = 0;
      std::size_t it_since_last_rho_update = 0;
      for (; sol.iterations <= settings.max_iterations;
           ++sol.iterations, ++it_since_last_rho_update)
      {
        // Fit the Anderson acceleration to compute accelerated x and y iterates
        if (sol.iterations > 1)
        {
          if (wk.anderson_history.capacity() > 0)
          {
            wk.anderson_history.push_back(wk.x, wk.z, wk.z - wk.z_previous);
          }

          if (
            wk.anderson_history.capacity() == 0 //
            || wk.anderson_history.size() < wk.anderson_history.capacity())
          {
            wk.x_anderson = wk.x;
            wk.z_anderson = wk.z;
          }
          else
          {
            wk.anderson_history.fit();
            wk.anderson_history.getAcceleratedIterates(wk.x_anderson, wk.z_anderson);
          }
        }

        // store previous iterates
        // note: when Anderson capacity is < 2, x_anderson_ = x_
        wk.x_previous = wk.x_anderson;
        wk.y_previous = wk.y;
        wk.z_previous = wk.z_anderson;
        sol.complementarity = Scalar(0);

        // y-update, using Anderson iterate.
        // If update is worse in terms of primal feas, it is rejected and the default
        // ADMM iterates are used to compute the y-update.
        {
          PINOCCHIO_TRACY_ZONE_SCOPED_N(
            "ADMMConstraintSolverTpl::solve - loop computeConeProjection");
          wk.tmp = wk.x_anderson - wk.z_anderson / (settings.tau * wk.rho);
          internal::computeConeProjection(constraint_models, constraint_datas, wk.tmp, wk.y);

          wk.anderson_primal_feasibility_vector = wk.x_anderson - wk.y;
          anderson_primal_feasibility =
            wk.anderson_primal_feasibility_vector.template lpNorm<Eigen::Infinity>();

          if (wk.anderson_history.capacity() > 1 && sol.iterations > 1)
          {
            if (
              anderson_primal_feasibility >= anderson_previous_primal_feasibility //
              && wk.anderson_history.size() == wk.anderson_history.capacity())
            {
              // Reject Anderson iterate, accept default ADMM iterate instead.
              // Reset Anderson acceleration.
              wk.x_previous = wk.x;
              wk.z_previous = wk.z;
              wk.tmp = wk.x_previous - wk.z_previous / (settings.tau * wk.rho);
              internal::computeConeProjection(constraint_models, constraint_datas, wk.tmp, wk.y);

              wk.anderson_history.clear();

              wk.anderson_primal_feasibility_vector = wk.x_previous - wk.y;
              anderson_primal_feasibility =
                wk.anderson_primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
            }
          }
        }
        anderson_previous_primal_feasibility = anderson_primal_feasibility;

        if (settings.solve_ncp)
        {
          // s-update
          internal::computeDeSaxeCorrection(
            constraint_models, constraint_datas, wk.z_previous, wk.desaxce);
        }

        // default (non-accelerated) x-update
        {
          PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMConstraintSolverTpl::solve - loop solveInPlace");
          wk.rhs =
            -(g + wk.desaxce - (wk.rho * settings.tau) * wk.y
              - (wk.mu_prox * settings.tau_prox) * wk.x_previous - wk.z_previous);
          wk.x = wk.rhs;
          G.solveInPlace(wk.x);
        }
        if (settings.stat_record)
        {
          G.applyOnTheRight(wk.x, wk.tmp);
          Scalar linear_system_residual = (wk.tmp - wk.rhs).template lpNorm<Eigen::Infinity>();
          stats.linear_system_residual.push_back(linear_system_residual);

          G.solveInPlace(wk.tmp);
          Scalar linear_system_consistency = (wk.tmp - wk.x).template lpNorm<Eigen::Infinity>();
          stats.linear_system_consistency.push_back(linear_system_consistency);
        }

        // default (non-accelerated) z-update
        wk.tmp = wk.z_previous - (settings.tau * wk.rho) * (wk.x - wk.y);
        wk.z.noalias() =
          settings.dual_momentum * wk.z_previous + (Scalar(1) - settings.dual_momentum) * wk.tmp;

        // check termination criteria
        wk.primal_feasibility_vector = wk.x - wk.y;

        {
          VectorXs & dx = wk.tmp;
          dx = wk.x - wk.x_previous;
          dx_norm = dx.template lpNorm<Eigen::Infinity>(); // check relative progress on x
          wk.dual_feasibility_vector = dx;
        }

        {
          VectorXs & dy = wk.tmp;
          dy = wk.y - wk.y_previous;
          dy_norm = dy.template lpNorm<Eigen::Infinity>(); // check relative progress on y
        }

        {
          VectorXs & dz = wk.tmp;
          dz = wk.z - wk.z_previous;
          dz_norm = dz.template lpNorm<Eigen::Infinity>(); // check relative progress on z
        }

        // compute primal/dual feasibility and complementarity
        // --> these are used to check convergence of the algo
        sol.primal_feasibility = wk.primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
        sol.dual_feasibility = wk.dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
        sol.dual_feasibility =
          math::max(wk.mu_prox * settings.tau_prox, wk.rho * settings.tau) * sol.dual_feasibility;
        internal::computeConicComplementarity(
          constraint_models, constraint_datas, wk.z, wk.y, sol.complementarity);

        if (settings.stat_record)
        {
          G.applyOnTheRight(wk.y, wk.rhs);
          wk.rhs += g - prox_value * wk.y;
          if (settings.solve_ncp)
          {
            internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.rhs, wk.tmp);
            wk.rhs += wk.tmp;
          }

          internal::computeDualConeProjection(constraint_models, constraint_datas, wk.rhs, wk.tmp);
          wk.rhs -= wk.tmp;

          Scalar dual_feasibility_ncp = wk.rhs.template lpNorm<Eigen::Infinity>();

          stats.primal_feasibility.push_back(sol.primal_feasibility);
          stats.dual_feasibility.push_back(sol.dual_feasibility);
          stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
          stats.complementarity.push_back(sol.complementarity);
          stats.rho.push_back(wk.rho);
          stats.mu_prox.push_back(wk.mu_prox);
          stats.anderson_size.push_back(wk.anderson_history.size());
        }

        // Checking stopping residual
        x_norm_inf = wk.x.template lpNorm<Eigen::Infinity>();
        y_norm_inf = wk.y.template lpNorm<Eigen::Infinity>();
        z_norm_inf = wk.z.template lpNorm<Eigen::Infinity>();
        // -- absolute check
        if (
          check_expression_if_real<Scalar, false>(
            sol.complementarity <= settings.tol_complementarity)
          && check_expression_if_real<Scalar, false>(
            sol.dual_feasibility
            <= settings.tol_feasibility
                 + settings.tol_rel_feasibility * math::max(g_norm_inf, z_norm_inf))
          && check_expression_if_real<Scalar, false>(
            sol.primal_feasibility
            <= settings.tol_feasibility
                 + settings.tol_rel_feasibility * math::max(x_norm_inf, y_norm_inf)))
        {
          abs_prec_reached = true;
        }
        else
        {
          abs_prec_reached = false;
        }

        // -- relative check
        if (
          check_expression_if_real<Scalar, false>(
            dx_norm <= settings.tol_rel_feasibility * math::max(x_norm_inf, x_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dy_norm <= settings.tol_rel_feasibility * math::max(y_norm_inf, y_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dz_norm <= settings.tol_rel_feasibility * math::max(z_norm_inf, z_previous_norm_inf)))
        {
          rel_prec_reached = true;
        }
        else
        {
          rel_prec_reached = false;
        }

        if (abs_prec_reached || rel_prec_reached)
          break;

        // update rho if needed
        if (
          sol.delassus_decomposition_update_count < settings.max_delassus_decomposition_updates
          && it_since_last_rho_update >= settings.rho_min_update_frequency)
        {
          // Apply rho according to the primal_dual_ratio
          Scalar new_rho = wk.rho;
          switch (settings.admm_update_rule)
          {
          case (ADMMUpdateRule::SPECTRAL):
            admm_update_rule_container.spectral_rule.eval(
              sol.primal_feasibility, sol.dual_feasibility, new_rho);
            break;
          case (ADMMUpdateRule::OSQP):
            admm_update_rule_container.osqp_rule.eval(
              sol.primal_feasibility, sol.dual_feasibility, new_rho);
            break;
          case (ADMMUpdateRule::LINEAR):
            admm_update_rule_container.linear_rule.eval(
              sol.primal_feasibility, sol.dual_feasibility, new_rho);
            break;
          case (ADMMUpdateRule::CONSTANT):
            break;
          }

          // clamp rho a second time
          new_rho = math::max(math::min(new_rho, rho_max), rho_min);

          // apply a momentum strategy on rho defined by:
          new_rho = std::pow(wk.rho, settings.rho_momentum)
                    * std::pow(new_rho, Scalar(1) - settings.rho_momentum);

          // clamp rho a second time in case the new values is outside the bounds
          new_rho = math::max(math::min(new_rho, rho_max), rho_min);

          bool update_delassus_factorization = false;
          if (new_rho == wk.rho)
          { // No change of rho, so need to redo a factorization
            update_delassus_factorization = false;
          }
          else if (
            new_rho >= settings.rho_update_ratio * wk.rho
            || wk.rho >= settings.rho_update_ratio * new_rho)
          { // sufficient change of the rho value
            wk.rho = new_rho;
            switch (settings.admm_proximal_rule)
            {
            case (ADMMProximalRule::MANUAL):
              // don't update the mu_prox
              break;
            case (ADMMProximalRule::AUTOMATIC):
              wk.mu_prox = wk.rho;
              break;
            }
            it_since_last_rho_update = 0;
            update_delassus_factorization = true;
          }

          // Account for potential update of rho
          if (update_delassus_factorization)
          {
            prox_value = settings.tau_prox * wk.mu_prox + settings.tau * wk.rho;
            if (old_prox_value != prox_value)
            {
              PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMConstraintSolverTpl::solve - loop updateDamping");
              wk.rhs = VectorXs::Constant(np, prox_value);
              G.updateDamping(wk.rhs);
              sol.delassus_decomposition_update_count++;
              old_prox_value = prox_value;
            }
          }
        }

        x_previous_norm_inf = x_norm_inf;
        y_previous_norm_inf = y_norm_inf;
        z_previous_norm_inf = z_norm_inf;

      } // end ADMM main for loop

      // Save values of spectral update rule
      if (settings.admm_update_rule == ADMMUpdateRule::SPECTRAL)
      {
        wk.spectral_rho_power = ADMMSpectralUpdateRule::computeRhoPower(L, m, wk.rho);
      }
    }

    PINOCCHIO_EIGEN_MALLOC_ALLOWED();

#ifdef PINOCCHIO_WITH_HPP_FCL
    if (settings.measure_timings)
    {
      timer.stop();
    }
#endif // PINOCCHIO_WITH_HPP_FCL

    if (settings.stat_record)
    {
      stats.iterations = sol.iterations;
      stats.delassus_decomposition_update_count = sol.delassus_decomposition_update_count;
    }

    sol.x = wk.x;
    sol.y = wk.y;
    sol.z = wk.z;
    sol.desaxce = wk.desaxce;
    sol.rho = wk.rho;
    sol.spectral_rho_power = wk.spectral_rho_power;
    sol.mu_prox = wk.mu_prox;
    sol.converged = abs_prec_reached || rel_prec_reached;

    return sol.converged;
  }

  template<typename _Scalar>
  template<
    typename DelassusDerived,
    typename VectorLike,
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  void ADMMConstraintSolverTpl<_Scalar>::retrievePrimalDualGuess(
    DelassusOperatorBase<DelassusDerived> & delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    ADMMSolverWorkspace & workspace,
    const ADMMSolverSettings & settings)
  {
    // for easier access
    ADMMSolverWorkspace & wk = workspace;
    DelassusDerived & G = delassus.derived();

    if (settings.primal_guess)
    {
      if (settings.dual_guess)
      {
        wk.z = settings.dual_guess.value();
        if (settings.solve_ncp)
        {
          // Add De Saxé shift
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.z, wk.desaxce);
          wk.z += wk.desaxce;
        }
        wk.x = settings.primal_guess.value();
        internal::computeConeProjection(constraint_models, constraint_datas, wk.x, wk.y);
      }
      else
      {
        // Warm-start dual variable using primal guess
        wk.x = settings.primal_guess.value();
        internal::computeConeProjection(constraint_models, constraint_datas, wk.x, wk.y);
        G.applyOnTheRight(wk.y, wk.z);
        wk.z.noalias() += g - wk.y.cwiseProduct(G.getDamping());
        if (settings.solve_ncp)
        {
          // Add De Saxé shift
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.z, wk.desaxce);
          wk.z += wk.desaxce;
        }
      }
    }
    else
    {
      if (settings.dual_guess)
      {
        // Warm-start primal variable using dual guess
        wk.z = settings.dual_guess.value();
        if (settings.solve_ncp)
        {
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.z, wk.desaxce);
          wk.z += wk.desaxce;
        }
        wk.x = wk.z - g - wk.desaxce;
        G.solveInPlace(wk.x);
        internal::computeConeProjection(constraint_models, constraint_datas, wk.x, wk.y);
        // wk.y.setZero();
      }
      else
      {
        wk.x.setZero();
        wk.y.setZero();
        wk.z = g;
        if (settings.solve_ncp)
        {
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.z, wk.desaxce);
          wk.z += wk.desaxce;
        }
      }
    }
  }

  //   template<typename _Scalar>
  //   template<
  //     typename DelassusDerived,
  //     typename VectorLike,
  //     template<typename T> class Holder,
  //     typename ConstraintModel,
  //     typename ConstraintModelAllocator>
  //   bool ADMMConstraintSolverTpl<_Scalar>::solve(
  //     DelassusOperatorBase<DelassusDerived> & _delassus,
  //     const Eigen::MatrixBase<VectorLike> & g,
  //     const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> &
  //     constraint_models, const Scalar dt, const boost::optional<RefConstVectorXs> preconditioner,
  //     const boost::optional<RefConstVectorXs> primal_guess,
  //     const boost::optional<RefConstVectorXs> dual_guess,
  //     const bool solve_ncp,
  //     const ADMMUpdateRule admm_update_rule,
  //     const bool stat_record)
  //
  //   {
  //     // Unused for now
  //     PINOCCHIO_UNUSED_VARIABLE(dual_guess);
  //
  //     typedef ADMMSpectralUpdateRuleTpl<Scalar> ADMMSpectralUpdateRule;
  //     typedef ADMMLinearUpdateRuleTpl<Scalar> ADMMLinearUpdateRule;
  //
  //     typedef ADMMUpdateRuleContainerTpl<Scalar> ADMMUpdateRuleContainer;
  //
  //     typedef DelassusOperatorPreconditionedTpl<DelassusDerived, DiagonalPreconditioner>
  //       DelassusOperatorPreconditioned;
  //     DelassusDerived & delassus = _delassus.derived();
  //
  //     const Scalar mu_R = delassus.getCompliance().minCoeff();
  //     PINOCCHIO_CHECK_INPUT_ARGUMENT(dt >= Scalar(0), "dt should be positive.");
  //     PINOCCHIO_CHECK_INPUT_ARGUMENT(tau <= Scalar(1) && tau > Scalar(0), "tau should lie in
  //     ]0,1]."); PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_prox >= 0, "mu_prox should be positive.");
  //     PINOCCHIO_CHECK_INPUT_ARGUMENT(mu_R >= Scalar(0), "R should be a positive vector.");
  //     PINOCCHIO_CHECK_ARGUMENT_SIZE(delassus.getCompliance().size(), problem_size);
  //
  //     // First, we initialize the primal and dual variables
  //     int it = 0;
  //     cholesky_update_count = 0;
  //
  //     Scalar complementarity, dx_bar_norm, dy_bar_norm, dz_bar_norm, //
  //       primal_feasibility, dual_feasibility;
  //
  //     if (stat_record)
  //     {
  //       stats.reserve(this->max_it);
  //       stats.reset();
  //     }
  //
  //     // Then, we get the time_scaling_acc_to_constraints T from the constraints to construct gs,
  //     // which is g time-scaled depending on the formulation of each constraint: gs = T^{-1} * g.
  //     The
  //     // idea is that if we formulate a given constraint at the position/velocity/acceleration
  //     level,
  //     // we want to measure constraint satisfaction for this constraint at the same
  //     // position/velocity/acceleration level.
  //     // However, to take admm steps, we work at the (force, acceleration) level for all
  //     constraints.
  //     // In short:
  //     // -> gs is used to perform optimization steps (we typically work on min_x x^TGx + gs^Tx).
  //     // -> time_scaling_acc_to_constraints is used to check for constraint satisfaction (we
  //     typically
  //     // want TGx + g = 0). It allows to go from accelerations to the units of each constraints.
  //     This
  //     // way, x and y are always forces expressed in N.
  //     // -> time_scaling_constraints_to_pos similarly allows to go from the units of the
  //     constraints
  //     // to positions in m. Warning: this constraints time-scaling has (a priori) nothing to do
  //     with
  //     // the pre-conditioner.
  //     internal::getTimeScalingFromAccelerationToConstraints(
  //       constraint_models, dt, time_scaling_acc_to_constraints);
  //     internal::getTimeScalingFromConstraintsToPosition(
  //       time_scaling_acc_to_constraints, dt, time_scaling_constraints_to_pos);
  //     gs = g.array() / time_scaling_acc_to_constraints.array();
  //     const Scalar g_pos_norm_inf =
  //       (g.cwiseProduct(time_scaling_constraints_to_pos)).template lpNorm<Eigen::Infinity>();
  //
  //     // Initialize De Saxé shift to 0
  //     // For the CCP, there is no shift
  //     // For the NCP, the shift will be initialized using z
  //     s_.setZero();
  //
  //     // Initial update of the variables
  //     // Init x
  //     if (primal_guess)
  //     {
  //       x_ = primal_guess.get();
  //       PINOCCHIO_CHECK_ARGUMENT_SIZE(x_.size(), problem_size);
  //     }
  //     else
  //     {
  //       x_.setZero();
  //     }
  //
  //     // Retrieve the pre-conditioner
  //     if (preconditioner)
  //     {
  //       PINOCCHIO_CHECK_ARGUMENT_SIZE(preconditioner_.rows(), problem_size);
  //       PINOCCHIO_CHECK_INPUT_ARGUMENT(
  //         preconditioner_.getDiagonal().minCoeff() > Scalar(0),
  //         "Preconditioner should be a strictly positive vector.");
  //       preconditioner_.setDiagonal(preconditioner.get());
  //     }
  //     else
  //     {
  //       preconditioner_.setDiagonal(VectorXs::Ones(problem_size));
  //     }
  //
  //     // Init y
  //     internal::computeConeProjection(constraint_models, x_, y_);
  //
  //     // Init z -> z_ = (G + R) * y_ + g
  //     delassus.applyOnTheRight(y_, z_);
  //     z_ += gs;
  //     z_ -= y_.cwiseProduct(delassus.getDamping());
  //     if (solve_ncp)
  //     {
  //       internal::computeDeSaxeCorrection(constraint_models, z_, s_);
  //       z_ += s_; // Add De Saxé shift
  //     }
  //
  //     /*
  //     // Computing the convergence criterion of the initial guess
  //     primal_feasibility = 0; // always feasible because y is projected
  //
  //     // complementarity of the initial guess
  //     // NB: complementarity is computed between a force y_ (in N) and z_ which unit is that of
  //     the
  //     // constraint formulation level.
  //     rhs =
  //       z_.array() * time_scaling_acc_to_constraints.array(); // back to constraint formulation
  //       level
  //     complementarity = internal::computeConicComplementarity(constraint_models, rhs, y_);
  //
  //     // dual feasibility is computed in "position" on the z_ variable (and not on z_bar_).
  //     dual_feasibility_vector = rhs;
  //     internal::computeDualConeProjection(constraint_models, rhs, rhs);
  //     dual_feasibility_vector -= rhs;
  //     dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
  //     const Scalar absolute_residual_warm_start = math::max(complementarity, dual_feasibility);
  //
  //     dual_feasibility_vector.array() *= time_scaling_constraints_to_pos.array();
  //     dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
  //     this->absolute_residual = math::max(complementarity, dual_feasibility);
  //
  //     // Checking if the initial guess is better than 0.
  //     // if instead of the x_ initial guess, x_ is set to 0, then z_ = g.
  //     // -> we check how much constraints violation is induced by using g as the dual variable.
  //     // note: here we work with g and not gs, because we check for constraints violation at the
  //     // formulation level of each constraints.
  //     const Scalar absolute_residual_zero_guess =
  //       computeZeroInitialGuessMaxConstraintViolation(constraint_models, g);
  //
  //     if (absolute_residual_zero_guess < absolute_residual_warm_start)
  //     { // If true, this means that the zero value initial guess leads a better feasibility in
  //     the
  //       // sense of the constraints satisfaction.
  //       // So we set the primal variables to the 0 initial guess and the dual variable to g.
  //       x_.setZero();
  //       y_.setZero();
  //       z_ = gs;
  //       if (solve_ncp)
  //       {
  //         {
  //           PINOCCHIO_TRACY_ZONE_SCOPED_N(
  //             "ADMMConstraintSolverTpl::solve - second computeDeSaxeCorrection");
  //           internal::computeDeSaxeCorrection(constraint_models, z_, s_);
  //         }
  //         z_ += s_; // Add De Saxé shift
  //       }
  //       rhs = z_.array() * time_scaling_acc_to_constraints.array();
  //       dual_feasibility_vector = rhs;
  //       internal::computeDualConeProjection(constraint_models, rhs, rhs);
  //       dual_feasibility_vector -= rhs; // Dual feasibility vector for the new null guess
  //       dual_feasibility_vector.array() *= time_scaling_constraints_to_pos.array();
  //       // We set the new convergence criterion
  //       this->absolute_residual = absolute_residual_zero_guess;
  //     }
  //     // We test convergence
  //     bool abs_prec_reached = this->absolute_residual < this->absolute_precision;
  //     */
  //
  //     bool abs_prec_reached = false;
  //
  //     x_.setZero();
  //     y_.setZero();
  //     z_ = gs;
  //     if (solve_ncp) {
  //       z_ += s_;
  //     }
  //
  //     if (!abs_prec_reached)
  //     { // the initial guess is not solution of the problem so we run the ADMM algorithm
  //       // Applying the preconditioner to work on a problem with a better scaling
  //       DelassusOperatorPreconditioned G_bar(_delassus, preconditioner_);
  //       rhs = VectorXs::Constant(this->problem_size, mu_prox);
  //       G_bar.updateDamping(rhs);      // G_bar =  P*(G+R)*P + mu_prox*Id
  //       scaleDualSolution(gs, g_bar_); // g_bar = P * gs
  //       scalePrimalSolution(x_, x_bar_);
  //       scalePrimalSolution(y_, y_bar_);
  //       scaleDualSolution(z_, z_bar_);
  //
  //       // Setup ADMM update rules:
  //       // Before running ADMM, we compute the largest and smallest eigenvalues of delassus in
  //       order
  //       // to be able to use a spectral update rule for the proximal parameter (rho)
  //       // TODO should we evaluate the eigenvalues of G or Gbar ?
  //       Scalar L, m, rho;
  //       ADMMUpdateRuleContainer admm_update_rule_container;
  //       switch (admm_update_rule)
  //       {
  //       case (ADMMUpdateRule::SPECTRAL): {
  //         if (this->problem_size > 1)
  //         {
  //           PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMConstraintSolverTpl::solve - lanczos");
  //           m = rhs.minCoeff();
  //           this->lanczos_decomposition.compute(G_bar);
  //           L = ::pinocchio::computeLargestEigenvalue(this->lanczos_decomposition.Ts(), 1e-8);
  // #ifndef NDEBUG
  //           const bool enforce_symmetry = true;
  //           Eigen::MatrixXd delassus = G_bar.matrix(enforce_symmetry);
  //           Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(delassus);
  //           Eigen::VectorXd eigvals = solver.eigenvalues();
  //           // Scalar true_m = eigvals.minCoeff();
  //           Scalar true_L = eigvals.maxCoeff();
  //           PINOCCHIO_UNUSED_VARIABLE(true_L);
  //             //          if (true_m > 0)
  //             //          {
  //             //            assert(
  //             //              math::fabs((true_m - m) / math::max(true_m, m)) < 0.01
  //             //              && "true_m and m are too far apart.");
  //             //          }
  //             // assert(
  //             //   math::fabs((true_L - L) / math::max(true_L, L)) < 0.01
  //             //   && "true_L and L are too far apart.");
  // #endif // NDEBUG
  //         }
  //         else
  //         {
  //           // TODO adapt this with Gbar
  //           typedef Eigen::Matrix<Scalar, 1, 1> Vector1;
  //           const Vector1 G = delassus * Vector1::Constant(1);
  //           m = L = G.coeff(0);
  //         }
  //         admm_update_rule_container.spectral_rule =
  //           ADMMSpectralUpdateRule(ratio_primal_dual, L, m, rho_power_factor);
  //         rho = ADMMSpectralUpdateRule::computeRho(L, m, rho_power);
  //         break;
  //       }
  //       case (ADMMUpdateRule::LINEAR):
  //         admm_update_rule_container.linear_rule =
  //           ADMMLinearUpdateRule(ratio_primal_dual, linear_update_rule_factor);
  //         rho = this->rho; // use the rho value stored in the solver.
  //         break;
  //       case (ADMMUpdateRule::CONSTANT):
  //         rho = this->rho; // use the rho value stored in the solver.
  //         break;
  //       }
  //
  //       // clamp the rho
  //       rho = math::max(1e-8, rho);
  //
  //       PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
  //
  //       // Update the cholesky decomposition
  //       Scalar prox_value = mu_prox + tau * rho;
  //       rhs = VectorXs::Constant(this->problem_size, prox_value);
  //       G_bar.updateDamping(rhs);
  //       Scalar old_prox_value = prox_value;
  //       cholesky_update_count = 1;
  //
  //       is_initialized = true;
  //
  //       // End of Initialization phase
  //       abs_prec_reached = false;
  //       bool rel_prec_reached = false;
  //
  //       Scalar x_bar_norm_inf = x_bar_.template lpNorm<Eigen::Infinity>();
  //       Scalar y_bar_norm_inf = y_bar_.template lpNorm<Eigen::Infinity>();
  //       Scalar z_bar_norm_inf = z_bar_.template lpNorm<Eigen::Infinity>();
  //       Scalar x_bar_previous_norm_inf = x_bar_norm_inf;
  //       Scalar y_bar_previous_norm_inf = y_bar_norm_inf;
  //       Scalar z_bar_previous_norm_inf = z_bar_norm_inf;
  //       it = 1;
  // #ifdef PINOCCHIO_WITH_HPP_FCL
  //       timer.start();
  // #endif // PINOCCHIO_WITH_HPP_FCL
  //       for (; it <= Base::max_it; ++it)
  //       {
  //
  //         x_bar_previous = x_bar_;
  //         y_bar_previous = y_bar_;
  //         z_bar_previous = z_bar_;
  //         complementarity = Scalar(0);
  //
  //         if (solve_ncp)
  //         {
  //           // s-update
  //           internal::computeDeSaxeCorrection(constraint_models, z_bar_, s_bar_);
  //         }
  //
  //         // x-update
  //         rhs = -(g_bar_ + s_bar_ - (rho * tau) * y_bar_ - mu_prox * x_bar_ - z_bar_);
  //         {
  //           PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMConstraintSolverTpl::solve - loop solveInPlace");
  //           G_bar.solveInPlace(rhs);
  //         }
  //         x_bar_ = rhs;
  //
  //         // y-update
  //         rhs -= z_bar_ / (tau * rho);
  //         {
  //           PINOCCHIO_TRACY_ZONE_SCOPED_N(
  //             "ADMMConstraintSolverTpl::solve - loop computeScaledConeProjection");
  //           internal::computeScaledConeProjection(
  //             constraint_models, rhs, preconditioner_.getDiagonal(), y_bar_);
  //         }
  //
  //         // z-update
  //         z_bar_ -= (tau * rho) * (x_bar_ - y_bar_);
  //
  //         // check termination criteria
  //         primal_feasibility_vector_bar = x_bar_ - y_bar_;
  //
  //         {
  //           VectorXs & dx_bar = rhs;
  //           dx_bar = x_bar_ - x_bar_previous;
  //           dx_bar_norm =
  //             dx_bar.template lpNorm<Eigen::Infinity>(); // check relative progress on x_bar
  //           dual_feasibility_vector_bar = mu_prox * dx_bar;
  //         }
  //
  //         {
  //           VectorXs & dy_bar = rhs;
  //           dy_bar = y_bar_ - y_bar_previous;
  //           dy_bar_norm =
  //             dy_bar.template lpNorm<Eigen::Infinity>(); // check relative progress on y_bar
  //           dual_feasibility_vector_bar += (tau * rho) * dy_bar;
  //         }
  //
  //         {
  //           VectorXs & dz_bar = rhs;
  //           dz_bar = z_bar_ - z_bar_previous;
  //           dz_bar_norm =
  //             dz_bar.template lpNorm<Eigen::Infinity>(); // check relative progress on z_bar
  //         }
  //
  //         // We unscale the quantities to work with stopping criterion from the original
  //         (unscaled)
  //         // problem
  //         unscalePrimalSolution(primal_feasibility_vector_bar, primal_feasibility_vector);
  //         primal_feasibility = primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
  //         unscaleDualSolution(dual_feasibility_vector_bar, dual_feasibility_vector);
  //         const Scalar dual_feasibility_admm =
  //           dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
  //         dual_feasibility_vector.array() *= time_scaling_acc_to_constraints.array();
  //         const Scalar dual_feasibility_constraint =
  //           dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
  //         dual_feasibility_vector.array() *= time_scaling_constraints_to_pos.array();
  //         dual_feasibility = dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
  //         unscalePrimalSolution(y_bar_, y_);
  //         unscaleDualSolution(z_bar_, z_);
  //         rhs = z_.array() * time_scaling_acc_to_constraints.array();
  //         complementarity = internal::computeConicComplementarity(constraint_models, rhs, y_);
  //
  //         if (stat_record)
  //         {
  //           VectorXs tmp(rhs);
  //           G_bar.applyOnTheRight(y_bar_, rhs);
  //           rhs += g_bar_ - prox_value * y_bar_;
  //           unscaleDualSolution(rhs, tmp);
  //           if (solve_ncp)
  //           {
  //             internal::computeDeSaxeCorrection(constraint_models, tmp, rhs);
  //             tmp += rhs;
  //           }
  //
  //           tmp.array() *=
  //             time_scaling_acc_to_constraints.array(); // back to constraint formulation level
  //           rhs = tmp;
  //           internal::computeDualConeProjection(constraint_models, rhs, rhs);
  //           tmp -= rhs;
  //
  //           Scalar dual_feasibility_ncp = tmp.template lpNorm<Eigen::Infinity>();
  //
  //           stats.primal_feasibility.push_back(primal_feasibility);
  //           stats.dual_feasibility.push_back(dual_feasibility);
  //           stats.dual_feasibility_admm.push_back(dual_feasibility_admm);
  //           stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
  //           stats.dual_feasibility_constraint.push_back(dual_feasibility_constraint);
  //           stats.complementarity.push_back(complementarity);
  //           stats.rho.push_back(rho);
  //         }
  //
  //         // Checking stopping residual
  //         const Scalar x_norm_inf = x_.template lpNorm<Eigen::Infinity>();
  //         const Scalar y_norm_inf = y_.template lpNorm<Eigen::Infinity>();
  //         const Scalar z_norm_inf = z_.template lpNorm<Eigen::Infinity>();
  //         if (
  //           check_expression_if_real<Scalar, false>(complementarity <= this->absolute_precision)
  //           && check_expression_if_real<Scalar, false>(
  //             dual_feasibility
  //             <= this->absolute_precision
  //                  + this->relative_precision * math::max(g_pos_norm_inf, z_norm_inf))
  //           && check_expression_if_real<Scalar, false>(
  //             primal_feasibility <= this->absolute_precision
  //                                     + this->relative_precision * math::max(x_norm_inf,
  //                                     y_norm_inf)))
  //           abs_prec_reached = true;
  //         else
  //           abs_prec_reached = false;
  //
  //         x_bar_norm_inf = x_bar_.template lpNorm<Eigen::Infinity>();
  //         y_bar_norm_inf = y_bar_.template lpNorm<Eigen::Infinity>();
  //         z_bar_norm_inf = z_bar_.template lpNorm<Eigen::Infinity>();
  //         if (
  //           check_expression_if_real<Scalar, false>(
  //             dx_bar_norm
  //             <= this->relative_precision * math::max(x_bar_norm_inf, x_bar_previous_norm_inf))
  //           && check_expression_if_real<Scalar, false>(
  //             dy_bar_norm
  //             <= this->relative_precision * math::max(y_bar_norm_inf, y_bar_previous_norm_inf))
  //           && check_expression_if_real<Scalar, false>(
  //             dz_bar_norm
  //             <= this->relative_precision * math::max(z_bar_norm_inf, z_bar_previous_norm_inf)))
  //           rel_prec_reached = true;
  //         else
  //           rel_prec_reached = false;
  //
  //         if (abs_prec_reached || rel_prec_reached)
  //           break;
  //
  //         // Apply rho according to the primal_dual_ratio
  //         bool update_delassus_factorization = false;
  //         switch (admm_update_rule)
  //         {
  //         case (ADMMUpdateRule::SPECTRAL):
  //           update_delassus_factorization = admm_update_rule_container.spectral_rule.eval(
  //             primal_feasibility, dual_feasibility_constraint, rho);
  //           break;
  //         case (ADMMUpdateRule::LINEAR):
  //           update_delassus_factorization = admm_update_rule_container.linear_rule.eval(
  //             primal_feasibility, dual_feasibility_constraint, rho);
  //           break;
  //         case (ADMMUpdateRule::CONSTANT):
  //           break;
  //         }
  //
  //         // clamp rho
  //         rho = math::max(1e-8, rho);
  //
  //         // Account for potential update of rho
  //         if (update_delassus_factorization)
  //         {
  //           prox_value = mu_prox + tau * rho;
  //           if (old_prox_value != prox_value)
  //           {
  //             PINOCCHIO_TRACY_ZONE_SCOPED_N("ADMMConstraintSolverTpl::solve - loop
  //             updateDamping"); rhs = VectorXs::Constant(this->problem_size, prox_value);
  //             G_bar.updateDamping(rhs);
  //             cholesky_update_count++;
  //             old_prox_value = prox_value;
  //           }
  //         }
  //
  //         x_bar_previous_norm_inf = x_bar_norm_inf;
  //         y_bar_previous_norm_inf = y_bar_norm_inf;
  //         z_bar_previous_norm_inf = z_bar_norm_inf;
  //       } // end ADMM main for loop
  //
  //       unscalePrimalSolution(x_bar_, x_);
  //       unscalePrimalSolution(y_bar_, y_);
  //       unscaleDualSolution(z_bar_, z_);
  //       unscaleDualSolution(s_bar_, s_);
  //
  //       this->relative_residual = math::max(
  //         dx_bar_norm / math::max(x_bar_norm_inf, x_bar_previous_norm_inf),
  //         dy_bar_norm / math::max(y_bar_norm_inf, y_bar_previous_norm_inf));
  //       this->relative_residual = math::max(
  //         this->relative_residual, dz_bar_norm / math::max(z_bar_norm_inf,
  //         z_bar_previous_norm_inf));
  //       this->absolute_residual =
  //         math::max(primal_feasibility, math::max(complementarity, dual_feasibility));
  //
  //       // Save values of spectral update rule
  //       if (admm_update_rule == ADMMUpdateRule::SPECTRAL)
  //       {
  //         this->rho_power = ADMMSpectralUpdateRule::computeRhoPower(L, m, rho);
  //         this->rho = rho;
  //       }
  //     }
  //     PINOCCHIO_EIGEN_MALLOC_ALLOWED();
  //
  // #ifdef PINOCCHIO_WITH_HPP_FCL
  //     timer.stop();
  // #endif // PINOCCHIO_WITH_HPP_FCL
  //     //
  //
  //     this->it = it;
  //     // we time-rescale dual solution and desaxe correction
  //     // so that z_ and s_ are back at the constraints formulations levels
  //     z_constraint_ = z_.array() * time_scaling_acc_to_constraints.array();
  //     s_constraint_ = s_.array() * time_scaling_acc_to_constraints.array();
  //
  //     if (stat_record)
  //     {
  //       stats.it = it;
  //       stats.cholesky_update_count = cholesky_update_count;
  //     }
  //
  //     if (abs_prec_reached)
  //       return true;
  //
  //     return false;
  //   }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_admm_solver_hxx__
