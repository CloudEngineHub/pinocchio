//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_solvers_clarabel_solver_hxx__
#define __pinocchio_algorithm_solvers_clarabel_solver_hxx__

#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT

  #include "pinocchio/macros.hpp"
  #include "pinocchio/algorithm/solvers/clarabel-solver.hpp"
  #include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"
  #include "pinocchio/utils/reference.hpp"

  #include <cpp/SupportedConeT.hpp>

namespace pinocchio
{

  namespace internal
  {

    template<typename VelocityVectorLike, typename ResultVectorLike>
    struct ClarabelDeSaxceShiftVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ClarabelDeSaxceShiftVisitor<VelocityVectorLike, ResultVectorLike>>
    {

      typedef boost::fusion::vector<const VelocityVectorLike &, ResultVectorLike &> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        ClarabelDeSaxceShiftVisitor<VelocityVectorLike, ResultVectorLike>>
        Base;
      using Base::run;

      template<typename ConstraintModel>
      static void algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const VelocityVectorLike & velocity,
        ResultVectorLike & result)
      {
        algo_impl(cmodel.derived(), velocity, result);
      }

      template<typename Scalar>
      static void algo_impl(
        const PointContactConstraintModelTpl<Scalar> & cmodel,
        const VelocityVectorLike & velocity,
        ResultVectorLike & result)
      {
        const Scalar mu = cmodel.getFriction();
        const Scalar vt_norm = velocity.template tail<2>().norm();
        // we follow clarabel convention
        result << mu * vt_norm, Scalar(0), Scalar(0);
      }

      template<typename ConstraintModel>
      static void algo_impl(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const VelocityVectorLike & velocity,
        ResultVectorLike & result)
      {
        PINOCCHIO_UNUSED_VARIABLE(cmodel);
        PINOCCHIO_UNUSED_VARIABLE(velocity);
        PINOCCHIO_UNUSED_VARIABLE(result);
        std::ostringstream msg;
        msg << cmodel.shortname() << " not yet supported in ClarabelConstraintSolver.";
        PINOCCHIO_THROW_PRETTY(std::runtime_error, msg.str());
      }

      template<typename ConstraintModel>
      static void run(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        const VelocityVectorLike & velocity,
        ResultVectorLike & result)
      {
        algo(cmodel.derived(), velocity, result);
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O> class ConstraintCollectionTpl>
      static void run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        const VelocityVectorLike & velocity,
        ResultVectorLike & result)
      {
        ArgsType args(velocity, result);
        run(cmodel.derived(), args);
      }
    }; // struct ClarabelDeSaxceShiftVisitor

    template<
      typename ConstraintModel,          //
      typename ConstraintModelAllocator, //
      typename ConstraintData,           //
      typename ConstraintDataAllocator,  //
      typename VectorLikeIn,             //
      typename VectorLikeOut>
    void computeClarabelDeSaxceShift(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const Eigen::MatrixBase<VectorLikeIn> & velocity,
      const Eigen::MatrixBase<VectorLikeOut> & shift_)
    {
      assert(velocity.size() == shift_.size());
      VectorLikeOut & shift = shift_.const_cast_derived();

      typedef typename VectorLikeIn::ConstSegmentReturnType SegmentType1;
      typedef typename VectorLikeOut::SegmentReturnType SegmentType2;

      Eigen::Index cindex = 0;
      for (std::size_t i = 0; i < constraint_models.size(); ++i)
      {
        const auto & cmodel = helper::get_ref(constraint_models[i]);
        const auto & cdata = helper::get_ref(constraint_datas[i]);
        const auto csize = cmodel.residualSize();

        SegmentType1 velocity_segment = velocity.segment(cindex, csize);
        SegmentType2 shift_segment = shift.segment(cindex, csize);

        typedef ClarabelDeSaxceShiftVisitor<SegmentType1, SegmentType2> Func;
        Func::run(cmodel, velocity_segment, shift_segment);

        cindex += csize;
      }
    }

    template<typename Scalar, typename FrictionVectorLike>
    struct ConstructClarabelConesVisitor
    : visitors::ConstraintUnaryVisitorBase<
        ConstructClarabelConesVisitor<Scalar, FrictionVectorLike>>
    {

      typedef std::vector<::clarabel::SupportedConeT<Scalar>> ClarabelConeVector;
      typedef boost::fusion::vector<ClarabelConeVector &, FrictionVectorLike &> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        ConstructClarabelConesVisitor<Scalar, FrictionVectorLike>>
        Base;
      using Base::run;

      template<typename ConstraintModel>
      static void algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        ClarabelConeVector & clarabel_cones,
        FrictionVectorLike & friction)
      {
        algo_impl(cmodel.derived(), clarabel_cones, friction);
      }

      static void algo_impl(
        const PointContactConstraintModelTpl<Scalar> & cmodel,
        ClarabelConeVector & clarabel_cones,
        FrictionVectorLike & friction)
      {
        const Scalar mu = cmodel.getFriction();
        // we follow clarabel convention
        friction << -mu, -Scalar(1), -Scalar(1);
        clarabel_cones.emplace_back(clarabel::SecondOrderConeT<Scalar>(3));
      }

      template<typename ConstraintModel>
      static void algo_impl(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        ClarabelConeVector & clarabel_cones,
        FrictionVectorLike & friction)
      {
        PINOCCHIO_UNUSED_VARIABLE(cmodel);
        PINOCCHIO_UNUSED_VARIABLE(clarabel_cones);
        PINOCCHIO_UNUSED_VARIABLE(friction);
        std::ostringstream msg;
        msg << cmodel.shortname() << " not yet supported in ClarabelConstraintSolver.";
        PINOCCHIO_THROW_PRETTY(std::runtime_error, msg.str());
      }

      template<typename ConstraintModel>
      static void run(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        ClarabelConeVector & clarabel_cones,
        FrictionVectorLike & friction)
      {
        algo(cmodel.derived(), clarabel_cones, friction);
      }

      template<int Options, template<typename S, int O> class ConstraintCollectionTpl>
      static void run(
        const ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        ClarabelConeVector & clarabel_cones,
        FrictionVectorLike & friction)
      {
        ArgsType args(clarabel_cones, friction);
        run(cmodel.derived(), args);
      }
    }; // struct ClarabelDeSaxceShiftVisitor

    template<
      typename Scalar,                   //
      typename ConstraintModel,          //
      typename ConstraintModelAllocator, //
      typename ConstraintData,           //
      typename ConstraintDataAllocator,  //
      typename VectorLikeFrictions>
    void constructClarabelCones(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      std::vector<::clarabel::SupportedConeT<Scalar>> & clarabel_cones,
      Eigen::MatrixBase<VectorLikeFrictions> & frictions_)
    {
      clarabel_cones.clear();
      clarabel_cones.reserve(constraint_models.size());
      VectorLikeFrictions & frictions = frictions_.const_cast_derived();

      typedef typename VectorLikeFrictions::SegmentReturnType SegmentType;

      Eigen::Index cindex = 0;
      for (std::size_t i = 0; i < constraint_models.size(); ++i)
      {
        const auto & cmodel = helper::get_ref(constraint_models[i]);
        const auto & cdata = helper::get_ref(constraint_datas[i]);
        const auto csize = cmodel.residualSize();

        SegmentType friction_segment = frictions.segment(cindex, csize);

        typedef ConstructClarabelConesVisitor<Scalar, SegmentType> Func;
        Func::run(cmodel, clarabel_cones, friction_segment);

        cindex += csize;
      }
    }

  } // namespace internal

  template<typename Scalar, int Options>
  template<
    typename DelassusDerived,
    typename VectorLike,
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  bool ClarabelConstraintSolverTpl<Scalar, Options>::solve(
    DelassusOperatorBase<DelassusDerived> & delassus,
    const Eigen::MatrixBase<VectorLike> & g,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const ClarabelSolverSettings & settings,
    ClarabelSolverResult & result)
  {
    // for easier access
    ClarabelSolverResult & res = result;
    ClarabelSolverWorkspace & ws = workspace_;

    {
      DelassusDerived & G = delassus.derived();
      static const Scalar max_init_numerical_damping = Eigen::NumTraits<Scalar>::dummy_precision();
      PINOCCHIO_THROW_IF(
        G.getDamping().minCoeff() < 0, std::runtime_error,
        "The delassus damping vector should be positive.");
      PINOCCHIO_THROW_IF(
        G.getDamping().maxCoeff() >= max_init_numerical_damping, std::runtime_error,
        "The delassus damping vector should be 0. Please call delassus.updateDamping(0) before "
        "calling the clarabel solver.");
    }

    const std::size_t problem_size = static_cast<std::size_t>(g.size());

    // -- check if settings are valid
    settings.checkValidity();

    // -- reset workspace
    // -- note: for now we only deal with point contact constraint so primal, dual and slack sizes
    // match
    ws.reset(problem_size);
    assert(g.size() == static_cast<Eigen::Index>(problem_size) && "Drift vector size mismatch");

    // Resize result vectors if needed
    res.reset(problem_size);

    // Build problem matrices for Clarabel
    // The contact problem is: min (1/2) lambda^T G lambda + g^T lambda
    //                         s.t. lambda in K (cone constraints)
    // In Simple, the convention is that the second-order cone is {(x, y, z), sqrt(x^2 + y^2) <= z}.
    // In Clarabel, the second-order cone convention is {(x, y, z), sqrt(y^2 + z^2) <= x}
    // So we have to do a permutation.
    // Get the dense Delassus matrix and the free velocity g
    using MapVectorXs = Eigen::Map<VectorXs>;
    using MapMatrixXs = Eigen::Map<MatrixXs>;
    const Eigen::Index np = static_cast<Eigen::Index>(problem_size);
    // we use alloca for a heap-allocation-free implementation
    MapMatrixXs G_dense(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, np, np));
    MapVectorXs q(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, np, 1));
    MapVectorXs q_new(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, np, 1));
    MapVectorXs b(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, np, 1));
    MapVectorXs frictions(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, np, 1));
    G_dense = delassus.derived().matrix(true); // true = enforce symmetry
    // permutate columns
    for (Eigen::Index i = 0; i < G_dense.cols(); i += 3)
    {
      G_dense.col(i).swap(G_dense.col(i + 1));
      G_dense.col(i).swap(G_dense.col(i + 2));
    }
    // permutate rows
    for (Eigen::Index i = 0; i < G_dense.rows(); i += 3)
    {
      G_dense.row(i).swap(G_dense.row(i + 1));
      G_dense.row(i).swap(G_dense.row(i + 2));
    }
    // Convert to sparse format for Clarabel
    // SparseMatrix P(G_dense.sparseView()); // = G_dense.sparseView();
    SparseMatrix P(G_dense.sparseView()); // = G_dense.sparseView();
    q = g;
    for (Eigen::Index i = 0; i < q.size(); i += 3)
    {
      std::swap(q.coeffRef(i), q.coeffRef(i + 1));
      std::swap(q.coeffRef(i), q.coeffRef(i + 2));
    }

    // Build cones and constraint matrix A and vector b
    // For the contact problem formulation, we have:
    // A * lambda + s = b, where s in K (friction cone)
    // We encode the friction coefficient in the diagonal A matrix
    b.setZero();
    internal::constructClarabelCones(
      constraint_models, constraint_datas, ws.clarabel_cones, frictions);
    // Apply permutation to friction coefficients to match Clarabel's cone convention
    SparseMatrix A(np, np);
    A.reserve(Eigen::VectorXi::Constant(np, 1));
    for (Eigen::Index i = 0; i < np; ++i)
    {
      const Scalar value = frictions(i);
      A.coeffRef(i, i) = value;
    }
    A.makeCompressed();

    // Create Clarabel settings from input
    clarabel::DefaultSettings<Scalar> clarabel_settings =
      clarabel::DefaultSettings<Scalar>::default_settings();
    clarabel_settings.max_iter = static_cast<uint32_t>(settings.max_iterations);
    clarabel_settings.tol_gap_abs = settings.absolute_complementarity_tol;
    clarabel_settings.tol_gap_rel = settings.relative_complementarity_tol;
    clarabel_settings.tol_feas = settings.absolute_feasibility_tol;
    clarabel_settings.tol_infeas_abs = settings.absolute_feasibility_tol;
    clarabel_settings.tol_infeas_rel = settings.relative_feasibility_tol;
    clarabel_settings.tol_ktratio = settings.tol_ktratio;
    clarabel_settings.time_limit = 100000;
    clarabel_settings.verbose = settings.verbose;

    // Create the Clarabel solver
    ws.clarabel_solver = std::make_shared<clarabel::DefaultSolver<Scalar>>(
      P, q, A, b, ws.clarabel_cones, clarabel_settings);

    // Solve the problem
    ws.clarabel_solver->solve();

    // Track iterations
    res.iterations = static_cast<std::size_t>(ws.clarabel_solver->info().iterations);

    ws.desaxce.setZero();
    if (settings.solve_ncp)
    {
      const clarabel::DefaultSolution<Scalar> & solution = ws.clarabel_solver->solution();
      ws.sigma.noalias() = G_dense * solution.x + q;
      internal::computeClarabelDeSaxceShift(
        constraint_models, constraint_datas, ws.sigma, ws.desaxce);
    }

    bool match = !settings.solve_ncp || (ws.desaxce.norm() <= settings.absolute_feasibility_tol);
    std::size_t ncp_loops = 1;
    ws.prev_desaxce = ws.desaxce;
    while (!match && ncp_loops <= settings.max_ncp_loops)
    {
      assert(settings.solve_ncp == true);
      q_new = q + ws.prev_desaxce;
      ws.clarabel_solver->update_q(q_new);
      ws.clarabel_solver->solve();
      res.iterations += static_cast<std::size_t>(ws.clarabel_solver->info().iterations);

      const clarabel::DefaultSolution<Scalar> & solution = ws.clarabel_solver->solution();
      ws.sigma.noalias() = G_dense * solution.x + q;
      internal::computeClarabelDeSaxceShift(
        constraint_models, constraint_datas, ws.sigma, ws.desaxce);
      if ((ws.desaxce - ws.prev_desaxce).norm() <= settings.absolute_feasibility_tol) // NCP measure
      {
        match = true;
      }
      ws.prev_desaxce = ws.desaxce;
      ++ncp_loops;
    }

    // Extract solution and store in result
    const clarabel::DefaultSolution<Scalar> & solution = ws.clarabel_solver->solution();

    // Store the solution with correct permutation for Simple's convention
    for (Eigen::Index i = 0; i < solution.x.size(); i += 3)
    {
      res.x(i) = solution.x(i + 1);
      res.x(i + 1) = solution.x(i + 2);
      res.x(i + 2) = solution.x(i);
      //
      res.z(i) = solution.z(i + 1);
      res.z(i + 1) = solution.z(i + 2);
      res.z(i + 2) = solution.z(i);
      //
      res.desaxce(i) = ws.desaxce(i + 1);
      res.desaxce(i + 1) = ws.desaxce(i + 2);
      res.desaxce(i + 2) = ws.desaxce(i);
      //
      res.sigma(i) = ws.sigma(i + 1);
      res.sigma(i + 1) = ws.sigma(i + 2);
      res.sigma(i + 2) = ws.sigma(i);
    }

    // Determine convergence based on Clarabel's solver status
    const bool has_converged =
      (ws.clarabel_solver->info().status == clarabel::SolverStatus::Solved);

    // Update result with convergence info
    res.converged = has_converged;
    res.primal_feasibility = solution.r_prim;
    res.dual_feasibility = solution.r_dual;
    res.complementarity = std::abs(res.x.dot(res.z));

    // Mark solver and result as valid
    is_valid_ = true;
    res.unsafe().makeValid();

    return res.converged;
  }

} // namespace pinocchio

#endif

#endif // __pinocchio_algorithm_solvers_clarabel_solver_hxx__
