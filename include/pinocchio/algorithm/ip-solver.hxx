//
// Copyright (c) 2022-2024 INRIA, KU Leuven
//

#ifndef __pinocchio_algorithm_ip_solver_hxx__
#define __pinocchio_algorithm_ip_solver_hxx__

#include "pinocchio/algorithm/contact-solver-utils.hpp"
#include "pinocchio/macros.hpp"
#include "pinocchio/utils/reference.hpp"

#include <iomanip>
#include <iostream>

#ifdef PINOCCHIO_WITH_HPP_FCL
  #define TIMER_START timer.start()
  #define TIMER_STOP timer.stop()
#else
  #define TIMER_START
  #define TIMER_STOP
#endif

namespace pinocchio
{
  // implementation
  template<typename Scalar>
  template<
    typename DelassusDerived,
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  void IPConstraintSolverTpl<Scalar>::solvePDSystem(
    const DelassusOperatorBase<DelassusDerived> & delassus,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    int iterative_refinement_steps)
  {
    iterative_refinement_steps = math::max(0, iterative_refinement_steps);
    rhs_x2 = rhs_x;
    rhs_z2 = rhs_z;
    rhs_s2 = rhs_s;
    delta_x.setZero();
    delta_s.setZero();
    delta_z.setZero();
    for (int k = 0; k <= iterative_refinement_steps; ++k)
    {
      ConeOps::inverseConeProduct(lambda, rhs_s, tmp_vec_0);
      tmp_vec_0 *= -1;
      tmp_vec_1 = rhs_z;
      for (Eigen::Index i = 0; i < problem_size / 3; ++i)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        tmp_vec_1.template segment<3>(3 * i).noalias() +=
          scaling.apply(tmp_vec_0.template segment<3>(3 * i));
      }
      // resulting gradient term addition (G^T W^{-1}W^{-T}g)
      // todo check if mem alloc
      // tmp_vec_2 is the barrier gradient
      auto & p_w_barrier = tmp_vec_2;
      auto & g = tmp_vec_1;
      for (Eigen::Index i = 0; i < problem_size / 3; ++i)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        Vector3s tmp = scaling.applyInverse(g.template segment<3>(3 * i));
        p_w_barrier.template segment<3>(3 * i).noalias() = -scaling.applyInverse(tmp);
      }
      denormalizeConeVariables(constraint_models, constraint_datas, p_w_barrier);
      p_w_barrier += rhs_x;
      delta_xi = -p_w_barrier;
      delassus.solveInPlace(delta_xi);
      // for z we have Wz = W^{-T})(g - G^T x)
      delta_zi.noalias() = -delta_xi;
      denormalizeConeVariables(constraint_models, constraint_datas, delta_zi);
      for (Eigen::Index i = 0; i < problem_size / 3; ++i)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        Vector3s tmp = g.template segment<3>(3 * i) + delta_zi.template segment<3>(3 * i);
        Vector3s tmp2 = scaling.applyInverse(tmp);
        delta_zi.template segment<3>(3 * i).noalias() = tmp2;
      }

      delta_si.noalias() = tmp_vec_0 - delta_zi;
      delta_x += delta_xi;
      delta_s += delta_si;
      delta_z += delta_zi;
      // compute residuals
      if (k == iterative_refinement_steps)
        break;
      rhs_x = rhs_x2;
      rhs_z = rhs_z2;
      rhs_s = rhs_s2;
      // compute rhs_x
      tmp_vec_0 = -delta_z;
      for (Eigen::Index i = 0; i < problem_size / 3; ++i)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        Vector3s tmp = tmp_vec_0.template segment<3>(3 * i);
        Vector3s tmp2 = scaling.applyInverse(tmp);
        tmp_vec_0.template segment<3>(3 * i).noalias() = tmp2;
      }
      denormalizeConeVariables(constraint_models, constraint_datas, tmp_vec_0);
      rhs_x += tmp_vec_0;
      tmp_vec_0 = delta_x;
      delassus.applyOnTheRight(tmp_vec_0, tmp_vec_1);
      // tmp_vec_1 += R.cwiseProduct(delta_x);
      rhs_x += tmp_vec_1;
      // compute rhs_s
      ConeOps::coneProduct(lambda, delta_s, tmp_vec_0);
      ConeOps::coneProduct(lambda, delta_z, tmp_vec_1);
      rhs_s += tmp_vec_0 + tmp_vec_1;
      // compute rhs_z
      tmp_vec_0 = -delta_x;
      denormalizeConeVariables(constraint_models, constraint_datas, tmp_vec_0);
      rhs_z += tmp_vec_0;
      for (Eigen::Index i = 0; i < problem_size / 3; ++i)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        Vector3s tmp = delta_s.template segment<3>(3 * i);
        Vector3s tmp2 = scaling.apply(tmp);
        rhs_z.template segment<3>(3 * i).noalias() += tmp2;
      }
    }
  }

  template<typename Scalar>
  template<
    typename DelassusDerived,
    typename VectorLikeConstraintDrift,
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  bool IPConstraintSolverTpl<Scalar>::solve(
    DelassusOperatorBase<DelassusDerived> & delassus,
    const Eigen::MatrixBase<VectorLikeConstraintDrift> & g,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const boost::optional<RefConstVectorXs> primal_guess,
    const boost::optional<RefConstVectorXs> dual_guess,
    bool solve_ncp,
    bool stat_record,
    int iterative_refinement_steps,
    Scalar target_mu,
    bool verbose)
  {
    // check if the delassus damping vector is a zero vector
    // TODO: should we set the damping to 0 and restore the damping at the end of solve?
    assert(delassus.derived().getDamping().isZero() && "The damping vector should be zero.");
    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();

    TIMER_START;
    if (stat_record)
    {
      stats.reserve(this->max_it);
      stats.reset();
    }

    v_constraint_evaluated = false;
    is_initialized = false;
    saxce_corr.setZero();

    // --------------
    // Initialization
    // --------------
    // Init x
    if (primal_guess)
    {
      x = primal_guess.get();
      if (dual_guess)
      {
        z = dual_guess.get();
      }
      else
      {
        // warmstart dual variable using primal guess
        internal::computeConeProjection(constraint_models, constraint_datas, x, rhs_x);
        delassus.applyOnTheRight(rhs_x, z);
        z += g;
      }
      if (solve_ncp)
      {
        internal::computeDeSaxeCorrection(constraint_models, constraint_datas, z, saxce_corr);
        z += saxce_corr;
      }
    }
    else
    {
      if (dual_guess)
      {
        // Warm-start primal variable using dual guess
        z = dual_guess.get();
        if (solve_ncp)
        {
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, z, saxce_corr);
          z += saxce_corr;
        }
        x = z - g - saxce_corr;
        delassus.updateDamping(1e-8);
        delassus.solveInPlace(x);
      }
      else
      {
        x.setZero();
        z = g;
        if (solve_ncp)
        {
          internal::computeDeSaxeCorrection(constraint_models, constraint_datas, z, saxce_corr);
          z += saxce_corr;
        }
      }
    }
    internal::computeDualConeProjection(constraint_models, constraint_datas, z, z);
    normalizeConeVariables(constraint_models, constraint_datas, z);
    internal::computeConeProjection(constraint_models, constraint_datas, x, s);
    denormalizeConeVariables(constraint_models, constraint_datas, s);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(x.size(), problem_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(s.size(), problem_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(z.size(), problem_size);

    // check if we already have a solution
    {
      computeCurrentConstraintVelocity(delassus, g);
      if (solve_ncp)
      {
        pinocchio::internal::computeDeSaxeCorrection(
          constraint_models, constraint_datas, v_constraint, saxce_corr);
      }
      computePrimalFeasibilityVector(constraint_models, constraint_datas);
      computePrimalOptimalityVector(constraint_models, constraint_datas);
      Scalar primal_feas = ccp_primal_feas.template lpNorm<Eigen::Infinity>();
      Scalar primal_opt = ccp_primal_opt.template lpNorm<Eigen::Infinity>()
                          / (std::max(1., (g + saxce_corr).template lpNorm<Eigen::Infinity>()));
      Scalar complementarity = s.dot(z);

      if (
        math::max(complementarity, math::max(primal_opt, primal_feas)) <= this->absolute_precision)
      {
        denormalizeConeVariables(constraint_models, constraint_datas, z);
        return true;
      }
    }

    // add a small value to the z components to make sure they are inside the
    // TODO: always initialize the solver, don't do the "rerun" thing
    // add a small value to the z components to make sure they are inside the
    // cone
    // Scalar min_dist = 1. * this->absolute_precision;
    Scalar min_dist = 1e-3;
    for (int i = 0; i < problem_size / 3; i++)
    {
      Scalar dist = s[3 * i + 2] - s.template segment<2>(3 * i).norm();
      if (dist < min_dist)
        s[3 * i + 2] = s.template segment<2>(3 * i).norm() + min_dist;
    }
    for (int i = 0; i < problem_size / 3; i++)
    {
      Scalar dist = z[3 * i + 2] - z.template segment<2>(3 * i).norm();
      if (dist < min_dist)
        z[3 * i + 2] = z.template segment<2>(3 * i).norm() + min_dist;
    }

    // initialize the lambda & scaling matrices
    for (int i = 0; i < problem_size / 3; i++)
    {
      auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
      lambda.template segment<3>(3 * i) =
        scaling.compute(s.template segment<3>(3 * i), z.template segment<3>(3 * i));
    }

    is_initialized = true;

    // ---------
    // Main loop
    // ---------
    if (verbose)
      printIterationsHeader();
    Scalar stepLength = 0.0;
    for (Base::it = 0; Base::it < Base::max_it; Base::it++)
    {
      //
      // Step 1 - Compute residuals
      //
      // constraint velocity
      computeCurrentConstraintVelocity(delassus, g);

      // de Saxcé correction term
      if (solve_ncp && Base::it % 5 == 0)
      {
        pinocchio::internal::computeDeSaxeCorrection(
          constraint_models, constraint_datas, v_constraint, saxce_corr);
      }

      // primal feasibility
      computePrimalFeasibilityVector(constraint_models, constraint_datas);

      // primal optimality
      computePrimalOptimalityVector(constraint_models, constraint_datas);

      // complementarity slackness
      computeComplementaritySlacknessVector();

      // (optional) objective value
      // delassus.applyOnTheRight(x, tmp_vec_0);
      // Scalar objective = 0.5 * x.dot(tmp_vec_0) + (g + saxce_corr).dot(x);
      // std::cout << "Objective value: " << objective << std::endl;

      //
      // Step 2 - Compute convergence precision
      //
      Scalar primal_feas = ccp_primal_feas.template lpNorm<Eigen::Infinity>();
      Scalar primal_opt = ccp_primal_opt.template lpNorm<Eigen::Infinity>()
                          / (std::max(1., (g + saxce_corr).template lpNorm<Eigen::Infinity>()));
      Scalar complementarity = 0.;
      for (int i = 0; i < problem_size / 3; i++)
      {
        complementarity =
          math::max(complementarity, math::fabs(ccp_compl_slackness[3 * i + 2] - target_mu));
      }

      if (verbose)
      {
        printIterationDetails(
          Base::it, primal_feas, primal_opt, complementarity, stepLength, complementarity);
      }

      if (stat_record)
      {
        stats.primal_feasibility.push_back(primal_feas);
        stats.dual_feasibility.push_back(primal_opt);
        // stats.dual_feasibility_ncp.push_back(dual_feasibility_ncp);
        stats.complementarity.push_back(complementarity);
        stats.mu.push_back(barrier_parameter);
      }

      //
      // Step 3 - Check if converged
      //
      if (
        (math::max(complementarity, primal_feas) <= this->absolute_precision)
        && primal_opt <= this->primal_opt_tol)
      {
        denormalizeConeVariables(constraint_models, constraint_datas, z);
        TIMER_STOP;
        if (stat_record)
        {
          stats.it = Base::it;
          stats.delassus_decomposition_update_count = Base::it;
        }
        PINOCCHIO_EIGEN_MALLOC_ALLOWED();
        return true; // TODO: replace with break
      }

      //
      // Step 4 - Compute search direction
      //
      // update the log-barrier parameter
      barrier_parameter = s.dot(z) / (problem_size / 3);

      // update the barrier terms
      // add the 3x3 block diagonal teems to the block diagonal of P matrix
      updateBarrierHessian(constraint_models, scaling_matrices, barrier_hessian_terms);

      // update delassus matrix and refactorize
      delassus.updateBarrierHessian(barrier_hessian_terms);

      // solve PD Newton systems
      // - first system  (i = 0): affine (predictor)
      // - second system (i = 1): corrected + centered (corrector + central path)
      Scalar dsdz;
      for (int i : {0, 1})
      {
        // set residuals
        if (i == 0)
        {
          if (barrier_parameter < target_mu)
          {
            continue;
          }
          rhs_x = ccp_primal_opt;
          rhs_z = ccp_primal_feas;
          rhs_s = ccp_compl_slackness;
          rhs_s(Eigen::seq(2, Eigen::indexing::last, 3)).array() -= target_mu;
        }
        else
        {
          rhs_x = ccp_primal_opt;
          rhs_z = ccp_primal_feas;
          rhs_s = ccp_compl_slackness;
          Scalar muplus;
          if (barrier_parameter > target_mu)
          {
            Scalar rho = 1. - stepLength + stepLength * stepLength * (dsdz / (lambda.dot(lambda)));
            Scalar sigma = math::pow(math::max(0., math::min(1., rho)), 3.);
            muplus = math::max(sigma * barrier_parameter, target_mu);
          }
          else
          {
            muplus = target_mu;
          }

          // add the perturbation to the rhs
          rhs_s(Eigen::seq(2, Eigen::indexing::last, 3)).array() -= muplus;

          // add the mehrotra correction to the rhs
          if (muplus > 10 * target_mu)
          {
            rhs_s += mehrotra_correction;
          }
        }

        // solve PD system for affine / combined direction
        solvePDSystem(delassus, constraint_models, constraint_datas, iterative_refinement_steps);
        if (i == 0)
        {
          // compute the mehrotra correction
          ConeOps::coneProduct(delta_s, delta_z, mehrotra_correction);

          // update dsdz
          dsdz = delta_s.dot(delta_z);
        }
        ConeOps::scale2(lambda, delta_s, delta_s);
        ConeOps::scale2(lambda, delta_z, delta_z);

        // compute the maximum step length stepLength
        Scalar t = 0.;
        for (int i = 0; i < problem_size / 3; i++)
        {
          t = math::max(t, delta_s.template segment<2>(3 * i).norm() - delta_s[3 * i + 2]);
          t = math::max(t, delta_z.template segment<2>(3 * i).norm() - delta_z[3 * i + 2]);
        }
        stepLength = (t == 0.) ? 1. : math::min(1., ((i == 0) ? 1.0 : 0.99) / t);
      }

      //
      // Step 6 - Take step
      //
      x = x + stepLength * delta_x;
      delta_s *= stepLength;
      delta_z *= stepLength;
      for (int i = 0; i < problem_size / 3; i++)
      {
        delta_s[3 * i + 2] += 1.;
        delta_z[3 * i + 2] += 1.;
      }
      ConeOps::scale2Inv(lambda, delta_s, delta_s);
      ConeOps::scale2Inv(lambda, delta_z, delta_z);
      for (int i = 0; i < problem_size / 3; i++)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        // s.template segment<3>(3*i) = scalings[i] * delta_s.template segment<3>(3*i);
        s.template segment<3>(3 * i) = scaling.apply(delta_s.template segment<3>(3 * i));
        z.template segment<3>(3 * i) = scaling.applyInverse(delta_z.template segment<3>(3 * i));
      }
      v_constraint_evaluated = false;

      //
      // Step 7 - Update lambda & scaling matrices -> move at the beginning of the loop
      //
      for (int i = 0; i < problem_size / 3; i++)
      {
        auto & scaling = scaling_matrices[static_cast<std::size_t>(i)];
        lambda.template segment<3>(3 * i) =
          scaling.update(delta_s.template segment<3>(3 * i), delta_z.template segment<3>(3 * i));
      }
    }
    TIMER_STOP
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    return false;
  };

  template<typename Scalar>
  void IPConstraintSolverTpl<Scalar>::printIterationDetails(
    int iteration,
    Scalar primal_feas,
    Scalar primal_opt,
    Scalar gap,
    Scalar step_length,
    Scalar ncp_convergence)
  {
    // Define column widths for formatting
    const int col_width = 20;

    // Reset fill to space for actual data
    std::cout << std::setfill(' ');

    // Print the values in the columns
    std::cout << std::setw(col_width) << iteration << std::setw(col_width) << std::scientific
              << primal_feas << std::setw(col_width) << std::scientific << primal_opt
              << std::setw(col_width) << std::scientific << gap << std::setw(col_width)
              << std::scientific << step_length << std::setw(col_width) << std::scientific
              << ncp_convergence << std::endl;
  }

  template<typename Scalar>
  void IPConstraintSolverTpl<Scalar>::printIterationsHeader()
  {
    // Define column widths for formatting
    const int col_width = 20;

    // Print the header
    std::cout << std::setw(col_width) << "Iteration" << std::setw(col_width) << "Primal Feasibility"
              << std::setw(col_width) << "Primal Optimality" << std::setw(col_width) << "Gap"
              << std::setw(col_width) << "Step Length" << std::setw(col_width) << "NCP Convergence"
              << std::endl;
  }

  template<typename VectorLikeInOut>
  struct NormalizeConeVariablesVisitor
  : visitors::ConstraintUnaryVisitorBase<NormalizeConeVariablesVisitor<VectorLikeInOut>>
  {
    typedef boost::fusion::vector<VectorLikeInOut &> ArgsType;
    typedef visitors::ConstraintUnaryVisitorBase<NormalizeConeVariablesVisitor<VectorLikeInOut>>
      Base;

    template<typename ConstraintModel>
    static void
    algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel, VectorLikeInOut & x)
    {
      PINOCCHIO_UNUSED_VARIABLE(cmodel);
      PINOCCHIO_UNUSED_VARIABLE(x);
      PINOCCHIO_THROW(
        std::runtime_error,
        "NormalizeConeVariablesVisitor not yet implemented for this constraint.");
    }

    template<typename Scalar, int Options>
    static void algo(
      const pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options> & cmodel,
      VectorLikeInOut & x)
    {
      x.coeffRef(2) /= cmodel.set().mu;
    }

    /// \brief Non-variant `run` method.
    template<typename ConstraintModel>
    static void
    run(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel, VectorLikeInOut & x)
    {
      algo(cmodel.derived(), x);
    }

    /// \brief Variant `run` method.
    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    static void run(
      const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      VectorLikeInOut & x)
    {
      ArgsType args(x);
      // Base::run will call `algo` (dispatch to right type)
      Base::run(cmodel, args);
    }
  };

  template<typename Scalar>
  template<
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator,
    typename VectorLikeInOut>
  void IPConstraintSolverTpl<Scalar>::normalizeConeVariables(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<VectorLikeInOut> & x_)
  {
    assert(
      constraint_models.size() == constraint_datas.size()
      && "Both std::vector should be of equal size.");

    using SegmentType = typename VectorLikeInOut::SegmentReturnType;

    VectorLikeInOut & x = x_.const_cast_derived();
    Eigen::DenseIndex cindex = 0;
    for (std::size_t k = 0; k < constraint_models.size(); ++k)
    {
      const auto & cmodel = helper::get_ref(constraint_models[k]);
      const auto & cdata = helper::get_ref(constraint_datas[k]);
      const auto csize = cmodel.activeSize(cdata);

      SegmentType x_segment = x.segment(cindex, csize);
      typedef NormalizeConeVariablesVisitor<SegmentType> Algo;
      Algo::run(cmodel, x_segment);

      cindex += csize;
    }
  }

  template<typename VectorLikeInOut>
  struct DenormalizeConeVariablesVisitor
  : visitors::ConstraintUnaryVisitorBase<DenormalizeConeVariablesVisitor<VectorLikeInOut>>
  {
    typedef boost::fusion::vector<VectorLikeInOut &> ArgsType;
    typedef visitors::ConstraintUnaryVisitorBase<DenormalizeConeVariablesVisitor<VectorLikeInOut>>
      Base;

    template<typename ConstraintModel>
    static void
    algo(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel, VectorLikeInOut & x)
    {
      PINOCCHIO_UNUSED_VARIABLE(cmodel);
      PINOCCHIO_UNUSED_VARIABLE(x);
      PINOCCHIO_THROW(
        std::runtime_error,
        "DenormalizeConeVariablesVisitor not yet implemented for this constraint.");
    }

    template<typename Scalar, int Options>
    static void algo(
      const pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options> & cmodel,
      VectorLikeInOut & x)
    {
      x.coeffRef(2) *= cmodel.set().mu;
    }

    /// \brief Non-variant `run` method.
    template<typename ConstraintModel>
    static void
    run(const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel, VectorLikeInOut & x)
    {
      algo(cmodel.derived(), x);
    }

    /// \brief Variant `run` method.
    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    static void run(
      const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      VectorLikeInOut & x)
    {
      ArgsType args(x);
      // Base::run will call `algo` (dispatch to right type)
      Base::run(cmodel, args);
    }
  };

  template<typename Scalar>
  template<
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator,
    typename VectorLikeInOut>
  void IPConstraintSolverTpl<Scalar>::denormalizeConeVariables(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<VectorLikeInOut> & x_)
  {
    assert(
      constraint_models.size() == constraint_datas.size()
      && "Both std::vector should be of equal size.");

    using SegmentType = typename VectorLikeInOut::SegmentReturnType;

    VectorLikeInOut & x = x_.const_cast_derived();
    Eigen::DenseIndex cindex = 0;
    for (std::size_t k = 0; k < constraint_models.size(); ++k)
    {
      const auto & cmodel = helper::get_ref(constraint_models[k]);
      const auto & cdata = helper::get_ref(constraint_datas[k]);
      const auto csize = cmodel.activeSize(cdata);

      SegmentType x_segment = x.segment(cindex, csize);
      typedef DenormalizeConeVariablesVisitor<SegmentType> Algo;
      Algo::run(cmodel, x_segment);

      cindex += csize;
    }
  }

  template<typename ScalingMatrix, typename BarrierHessianTerm>
  struct UpdateBarrierHessianTermVisitor
  : visitors::ConstraintUnaryVisitorBase<
      UpdateBarrierHessianTermVisitor<ScalingMatrix, BarrierHessianTerm>>
  {
    typedef boost::fusion::vector<const ScalingMatrix &, BarrierHessianTerm &> ArgsType;
    typedef visitors::ConstraintUnaryVisitorBase<
      UpdateBarrierHessianTermVisitor<ScalingMatrix, BarrierHessianTerm>>
      Base;

    template<typename ConstraintModel>
    static void algo(
      const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
      const ScalingMatrix & scaling_matrix,
      BarrierHessianTerm & barrier_hessian_term)
    {
      PINOCCHIO_UNUSED_VARIABLE(cmodel);
      PINOCCHIO_UNUSED_VARIABLE(scaling_matrix);
      PINOCCHIO_UNUSED_VARIABLE(barrier_hessian_term);
      PINOCCHIO_THROW(
        std::runtime_error,
        "UpdateBarrierHessianTermVisitor not yet implemented for this constraint.");
    }

    template<typename Scalar, int Options>
    static void algo(
      const pinocchio::FrictionalPointConstraintModelTpl<Scalar, Options> & cmodel,
      const ScalingMatrix & scaling_matrix,
      BarrierHessianTerm & barrier_hessian_term)
    {
      auto W_minTG = scaling_matrix.getInverseMatrix();
      W_minTG.col(2) *= cmodel.set().mu;
      barrier_hessian_term = W_minTG.transpose() * W_minTG;
    }

    /// \brief Non-variant `run` method.
    template<typename ConstraintModel>
    static void run(
      const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
      const ScalingMatrix & scaling_matrix,
      BarrierHessianTerm & barrier_hessian_term)
    {
      algo(cmodel.derived(), scaling_matrix, barrier_hessian_term);
    }

    /// \brief Variant `run` method.
    template<
      typename Scalar,
      int Options,
      template<typename S, int O> class ConstraintCollectionTpl>
    static void run(
      const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const ScalingMatrix & scaling_matrix,
      BarrierHessianTerm & barrier_hessian_term)
    {
      ArgsType args(scaling_matrix, barrier_hessian_term);
      // Base::run will call `algo` (dispatch to right type)
      Base::run(cmodel, args);
    }
  };

  template<typename Scalar>
  template<typename ConstraintModel, typename ConstraintModelAllocator>
  void IPConstraintSolverTpl<Scalar>::updateBarrierHessian(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const ScalingMatrixVector & scaling_matrices,
    BarrierHessianTermVector & barrier_hessian_terms)
  {
    assert(
      constraint_models.size() == scaling_matrices.size()
      && "Both std::vector should be of equal size.");
    assert(
      constraint_models.size() == barrier_hessian_terms.size()
      && "Both std::vector should be of equal size.");

    for (std::size_t k = 0; k < constraint_models.size(); ++k)
    {
      const auto & cmodel = helper::get_ref(constraint_models[k]);
      const auto & scaling_matrix = helper::get_ref(scaling_matrices[k]);
      auto & barrier_hessian_term = helper::get_ref(barrier_hessian_terms[k]);

      typedef UpdateBarrierHessianTermVisitor<ScalingMatrix, BarrierHessianTerm> Algo;
      Algo::run(cmodel, scaling_matrix, barrier_hessian_term);
    }
  }

  template<typename Scalar>
  template<
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  void IPConstraintSolverTpl<Scalar>::computePrimalFeasibilityVector(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    ccp_primal_feas.noalias() = -x;
    denormalizeConeVariables(constraint_models, constraint_datas, ccp_primal_feas);
    ccp_primal_feas.noalias() += s;
  }

  template<typename Scalar>
  template<
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  void IPConstraintSolverTpl<Scalar>::computePrimalOptimalityVector(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    ccp_primal_opt.noalias() = -z;
    denormalizeConeVariables(constraint_models, constraint_datas, ccp_primal_opt);
    ccp_primal_opt.noalias() += v_constraint;
    ccp_primal_opt.noalias() += saxce_corr;
  }

  template<typename Scalar>
  void IPConstraintSolverTpl<Scalar>::computeComplementaritySlacknessVector()
  {
    ConeOps::coneProduct(lambda, lambda, ccp_compl_slackness);
  }

} // namespace pinocchio

#undef TIMER_START
#undef TIMER_STOP

#endif // ifndef __pinocchio_algorithm_ip_solver_hpp__
