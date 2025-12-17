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
    const ADMMSolverSettings & settings,
    ADMMSolverResult & result)
  {
    // for easier access
    static constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();
    ADMMSolverResult & res = result;
    ADMMSolverWorkspace & wk = workspace_;
    DelassusDerived & G = delassus.derived();

    // Configure/reset workspace, stats and result
    // note: the order matters as workspace is initialized using
    // optional warmstarts contained in result.
    const Eigen::Index np = G.rows();
    const std::size_t problem_size = static_cast<std::size_t>(np);
    assert(G.cols() == np);
    assert(g.size() == np);
    assert(residualSize(constraint_models, constraint_datas) == np);

    // -- check if settings are valid
    settings.checkValidity();

    // -- reset workspace
    wk.reset(problem_size, settings.lanczos_size, settings.anderson_capacity);
    assert(wk.problem_size == problem_size);
    assert(wk.x.size() == np);

    // -- reset per-iteration statistics
    stats.reset();
    if (settings.stat_record)
    {
      stats.reserve(settings.max_iterations);
    }

    // -- retrieve warmstart from results, then reset results
    retrievePrimalDualGuess(delassus, g, constraint_models, constraint_datas, settings, res, wk);
    retrieveRhoGuess(delassus, settings, res, wk);
    switch (settings.admm_proximal_rule)
    {
    case (ADMMProximalRule::MANUAL):
      wk.mu_prox = settings.mu_prox;
      break;
    case (ADMMProximalRule::AUTOMATIC):
      wk.mu_prox = wk.rho;
      break;
    }
    res.reset(problem_size);
    assert(res.isValid() == false);
    assert(res.problem_size == problem_size);
    assert(res.iterations == 0);

    // -- init of internals done - the solver is now marked as reset
    is_valid_ = false;

#ifdef PINOCCHIO_WITH_HPP_FCL
    if (settings.measure_timings)
    {
      timer.start();
    }
#endif // PINOCCHIO_WITH_HPP_FCL

    // Check NCP/CCP conditions. If they are satisfied, don't run the solver.
    // -- always primaly feasible as y is projected onto the constraints
    // (this has been done in `retrievePrimalDualGuess`)
    res.primal_feasibility = Scalar(0);

    // -- dual feasibility
    G.applyOnTheRight(wk.y, wk.rhs);
    wk.rhs += g - wk.y.cwiseProduct(G.getDamping());
    if (settings.solve_ncp)
    {
      internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.rhs, wk.desaxce);
      wk.rhs += wk.desaxce;
    }
    else
    {
      wk.desaxce.setZero();
    }
    internal::computeDualConeProjection(constraint_models, constraint_datas, wk.rhs, wk.tmp);
    wk.tmp -= wk.rhs;
    res.dual_feasibility = wk.tmp.template lpNorm<Eigen::Infinity>();

    // -- complementarity
    internal::computeConicComplementarity(
      constraint_models, constraint_datas, wk.rhs, wk.y, res.complementarity);

    bool abs_prec_reached = false;
    bool rel_prec_reached = false;
    if (
      check_expression_if_real<Scalar, false>(
        res.complementarity <= settings.absolute_tol_complementarity)
      && check_expression_if_real<Scalar, false>(
        res.dual_feasibility <= settings.absolute_tol_feasibility))
    {
      abs_prec_reached = true;
      wk.z = wk.rhs; // store dual solution
    }

    if (!abs_prec_reached)
    {
      // Setup ADMM update rules:
      // Before running ADMM, we compute the largest and smallest eigenvalues of delassus in order
      // to be able to use a spectral update rule for the proximal parameter (rho)
      ADMMUpdateRuleContainer admm_update_rule_container;
      switch (settings.admm_update_rule)
      {
      case (ADMMUpdateRule::SPECTRAL): {
        if (wk.delassus_largest_eigenvalue.has_value() == false)
        {
          // largest eigenvalue has not yet been computed, we compute it
          wk.delassus_largest_eigenvalue = computeDelassusLargestEigenvalue(G, wk);
        }
        admm_update_rule_container.spectral_rule = ADMMSpectralUpdateRule(
          settings.ratio_primal_dual,              //
          wk.delassus_largest_eigenvalue.value(),  //
          wk.delassus_smallest_eigenvalue.value(), //
          settings.spectral_rho_power_factor);
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
      wk.rhs.setConstant(prox_value);
      G.updateDamping(wk.rhs);
      Scalar old_prox_value = prox_value;
      wk.delassus_decomposition_update_count++;

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

      res.iterations = 0;
      std::size_t it_since_last_rho_update = 0;
      for (; res.iterations <= settings.max_iterations;
           ++res.iterations, ++it_since_last_rho_update)
      {
        // Fit the Anderson acceleration to compute accelerated x and y iterates
        if (res.iterations > 1)
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
        res.complementarity = Scalar(0);

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

          if (wk.anderson_history.capacity() > 1 && res.iterations > 1)
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
        else
        {
          wk.desaxce.setZero();
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
          auto & dx = wk.tmp;
          dx = wk.x - wk.x_previous;
          dx_norm = dx.template lpNorm<Eigen::Infinity>(); // check relative progress on x
          wk.dual_feasibility_vector = dx;
        }

        {
          auto & dy = wk.tmp;
          dy = wk.y - wk.y_previous;
          dy_norm = dy.template lpNorm<Eigen::Infinity>(); // check relative progress on y
        }

        {
          auto & dz = wk.tmp;
          dz = wk.z - wk.z_previous;
          dz_norm = dz.template lpNorm<Eigen::Infinity>(); // check relative progress on z
        }

        // compute primal/dual feasibility and complementarity
        // --> these are used to check convergence of the algo
        res.primal_feasibility = wk.primal_feasibility_vector.template lpNorm<Eigen::Infinity>();
        res.dual_feasibility = wk.dual_feasibility_vector.template lpNorm<Eigen::Infinity>();
        res.dual_feasibility =
          math::max(wk.mu_prox * settings.tau_prox, wk.rho * settings.tau) * res.dual_feasibility;
        internal::computeConicComplementarity(
          constraint_models, constraint_datas, wk.z, wk.y, res.complementarity);

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

          stats.primal_feasibility.push_back(res.primal_feasibility);
          stats.dual_feasibility.push_back(res.dual_feasibility);
          stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
          stats.complementarity.push_back(res.complementarity);
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
            res.complementarity <= settings.absolute_tol_complementarity)
          && check_expression_if_real<Scalar, false>(
            res.dual_feasibility
            <= settings.absolute_tol_feasibility
                 + settings.relative_tol_feasibility * math::max(g_norm_inf, z_norm_inf))
          && check_expression_if_real<Scalar, false>(
            res.primal_feasibility
            <= settings.absolute_tol_feasibility
                 + settings.relative_tol_feasibility * math::max(x_norm_inf, y_norm_inf)))
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
            dx_norm
            <= settings.relative_tol_feasibility * math::max(x_norm_inf, x_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dy_norm
            <= settings.relative_tol_feasibility * math::max(y_norm_inf, y_previous_norm_inf))
          && check_expression_if_real<Scalar, false>(
            dz_norm
            <= settings.relative_tol_feasibility * math::max(z_norm_inf, z_previous_norm_inf)))
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
          wk.delassus_decomposition_update_count < settings.max_delassus_decomposition_updates
          && it_since_last_rho_update >= settings.rho_min_update_frequency)
        {
          // Apply rho according to the primal_dual_ratio
          Scalar new_rho = wk.rho;
          switch (settings.admm_update_rule)
          {
          case (ADMMUpdateRule::SPECTRAL):
            admm_update_rule_container.spectral_rule.eval(
              res.primal_feasibility, res.dual_feasibility, new_rho);
            break;
          case (ADMMUpdateRule::OSQP):
            admm_update_rule_container.osqp_rule.eval(
              res.primal_feasibility, res.dual_feasibility, new_rho);
            break;
          case (ADMMUpdateRule::LINEAR):
            admm_update_rule_container.linear_rule.eval(
              res.primal_feasibility, res.dual_feasibility, new_rho);
            break;
          case (ADMMUpdateRule::CONSTANT):
            break;
          }

          // clamp rho a second time
          new_rho = math::max(math::min(new_rho, settings.rho_max), settings.rho_min);

          // apply a momentum strategy on rho defined by:
          new_rho = std::pow(wk.rho, settings.rho_momentum)
                    * std::pow(new_rho, Scalar(1) - settings.rho_momentum);

          // clamp rho a second time in case the new values is outside the bounds
          new_rho = math::max(math::min(new_rho, settings.rho_max), settings.rho_min);

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
              wk.rhs.setConstant(prox_value);
              G.updateDamping(wk.rhs);
              wk.delassus_decomposition_update_count++;
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
        wk.spectral_rho_power = ADMMSpectralUpdateRule::computeRhoPower(
          wk.delassus_largest_eigenvalue.value(),  //
          wk.delassus_smallest_eigenvalue.value(), //
          wk.rho);
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
      stats.iterations = res.iterations;
      stats.delassus_decomposition_update_count = wk.delassus_decomposition_update_count;
    }

    res.x = wk.x;
    res.y = wk.y;
    res.z = wk.z;
    res.desaxce = wk.desaxce;
    res.rho = wk.rho;
    res.spectral_rho_power = wk.spectral_rho_power;
    res.mu_prox = wk.mu_prox;
    res.delassus_decomposition_update_count = wk.delassus_decomposition_update_count;
    res.converged = abs_prec_reached || rel_prec_reached;
    res.makeValid();

    // the solver has run, we mark it as valid
    is_valid_ = true;

    return res.converged;
  }

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
  void ADMMConstraintSolverTpl<_Scalar>::retrievePrimalDualGuess(
    DelassusOperatorBase<DelassusDerived> & delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const ADMMSolverSettings & settings,
    const ADMMSolverResult & result,
    ADMMSolverWorkspace & workspace)
  {
    // for easier access
    const ADMMSolverResult & res = result;
    DelassusDerived & G = delassus.derived();
    ADMMSolverWorkspace & wk = workspace;

    const Scalar min_compliance = G.getCompliance().minCoeff();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      min_compliance >= Scalar(0), "compliance should be a positive vector.");

    // initialize De Saxé shift to 0
    // for the CCP, there is no shift.
    // for the NCP, the shift will be initialized using z.
    wk.desaxce.setZero();

    // set initial damping of the delassus to the proximal value and get smallest possible
    // eigenvalue of the problem.
    wk.mu_prox = settings.mu_prox;
    wk.delassus_smallest_eigenvalue = min_compliance + wk.mu_prox;
    wk.rhs.setConstant(wk.mu_prox);
    G.updateDamping(wk.rhs);
    wk.delassus_decomposition_update_count++;

    // Initialization of the primal/dual variables.
    // If both primal and dual guesses are given, the solver uses both.
    // If one is given but not the other, the solver will compute the missing one using the given
    // one.
    if (res.primal_guess)
    {
      if (res.dual_guess)
      {
        wk.z = res.dual_guess.value();
        if (settings.solve_ncp)
        {
          // Add De Saxé shift
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, wk.z, wk.desaxce);
          wk.z += wk.desaxce;
        }
        wk.x = res.primal_guess.value();
        internal::computeConeProjection(constraint_models, constraint_datas, wk.x, wk.y);
      }
      else
      {
        // Warm-start dual variable using primal guess
        wk.x = res.primal_guess.value();
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
      if (res.dual_guess)
      {
        // Warm-start primal variable using dual guess
        wk.z = res.dual_guess.value();
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

    // sanity checks
    const Eigen::Index np = G.rows();
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.x.size(), np);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.y.size(), np);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.z.size(), np);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(wk.desaxce.size(), np);
  }

  template<typename _Scalar>
  template<typename DelassusDerived>
  void ADMMConstraintSolverTpl<_Scalar>::retrieveRhoGuess(
    const DelassusOperatorBase<DelassusDerived> & delassus,
    const ADMMSolverSettings & settings,
    const ADMMSolverResult & result,
    ADMMSolverWorkspace & workspace)
  {
    // for easier access
    const DelassusDerived & G = delassus.derived();
    const ADMMSolverResult & res = result;
    ADMMSolverWorkspace & wk = workspace;

    std::optional<Scalar> rho_init = settings.rho_init;
    Scalar spectral_rho_power_init = settings.spectral_rho_power_init;
    if (settings.warmstart_rho_with_prev_sol && res.isValid())
    {
      // override rho_init with previous result's rho value
      rho_init = res.rho;
      spectral_rho_power_init = res.spectral_rho_power;
    }

    // init workspace's rho parameters
    wk.spectral_rho_power = spectral_rho_power_init;
    if (rho_init)
    {
      wk.rho = rho_init.value();
    }
    else
    {
      // compute rho with spectral rule
      assert(wk.delassus_smallest_eigenvalue.has_value() == true);
      assert(wk.delassus_largest_eigenvalue.has_value() == false);
      wk.delassus_largest_eigenvalue = computeDelassusLargestEigenvalue(G, wk);
      if (std::isnan(wk.delassus_largest_eigenvalue.value()))
      {
        wk.delassus_largest_eigenvalue = Scalar(1);
      }
      // TODO: change order of largest/smallest
      wk.rho = ADMMSpectralUpdateRule::computeRho(
        wk.delassus_largest_eigenvalue.value(),  //
        wk.delassus_smallest_eigenvalue.value(), //
        wk.spectral_rho_power);
    }
    PINOCCHIO_CHECK_INPUT_ARGUMENT(wk.rho >= 0, "rho should be positive.");

    // clamp the rho
    wk.rho = math::max(math::min(wk.rho, settings.rho_max), settings.rho_min);
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_solvers_admm_solver_hxx__
