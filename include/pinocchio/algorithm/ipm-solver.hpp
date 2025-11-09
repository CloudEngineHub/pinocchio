//
// Copyright (c) 2022-2024 INRIA, KU Leuven
//

#ifndef __pinocchio_algorithm_ip_solver_hpp__
#define __pinocchio_algorithm_ip_solver_hpp__

#include "pinocchio/algorithm/ipm-solver-cone-operations.hpp"
#include "pinocchio/algorithm/contact-solver-base.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"

#include <boost/optional.hpp>

namespace pinocchio
{
  /// \brief Interior point solver
  template<typename _Scalar>
  struct IPMConstraintSolverTpl : ContactSolverBaseTpl<_Scalar>
  {
    typedef _Scalar Scalar;
    typedef ContactSolverBaseTpl<Scalar> Base;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3s;
    typedef Eigen::Matrix<Scalar, 3, 3> Matrix3x3;
    typedef CoulombFrictionConeTpl<Scalar> Cone;
    typedef IPMSolverConeOperations<Scalar> ConeOps;
    typedef Eigen::Ref<const VectorXs> RefConstVectorXs;
    typedef Matrix3x3 BarrierHessianTerm;
    typedef std::vector<BarrierHessianTerm> BarrierHessianTermVector;
    typedef typename ConeOps::ScalingMatrix ScalingMatrix;
    typedef std::vector<ScalingMatrix> ScalingMatrixVector;

    explicit IPMConstraintSolverTpl(const int problem_size)
    : Base(problem_size)
    , v_constraint(VectorXs::Zero(problem_size))
    , x(VectorXs::Zero(problem_size))
    , x_previous(VectorXs::Zero(problem_size))
    , z(VectorXs::Zero(problem_size))
    , s(VectorXs::Zero(problem_size))
    , lambda(VectorXs::Zero(problem_size))
    , delta_x(VectorXs::Zero(problem_size))
    , delta_z(VectorXs::Zero(problem_size))
    , delta_s(VectorXs::Zero(problem_size))
    , delta_xi(VectorXs::Zero(problem_size))
    , delta_zi(VectorXs::Zero(problem_size))
    , delta_si(VectorXs::Zero(problem_size))
    , ccp_primal_feas(VectorXs::Zero(problem_size))
    , ccp_primal_opt(VectorXs::Zero(problem_size))
    , ccp_compl_slackness(VectorXs::Zero(problem_size))
    , ncp_primal_opt(VectorXs::Zero(problem_size))
    , tmp_vec_0(VectorXs::Zero(problem_size))
    , tmp_vec_1(VectorXs::Zero(problem_size))
    , tmp_vec_2(VectorXs::Zero(problem_size))
    , mehrotra_correction(VectorXs::Zero(problem_size))
    , rhs_x(VectorXs::Zero(problem_size))
    , rhs_z(VectorXs::Zero(problem_size))
    , rhs_s(VectorXs::Zero(problem_size))
    , rhs_x2(VectorXs::Zero(problem_size))
    , rhs_z2(VectorXs::Zero(problem_size))
    , rhs_s2(VectorXs::Zero(problem_size))
    , saxce_corr(VectorXs::Zero(problem_size))
    , barrier_hessian_terms(std::size_t(problem_size / 3), Matrix3x3::Zero())
    , scaling_matrices(std::size_t(problem_size / 3))
    , stats(Base::max_it)
    {
    }

    struct IPMSolverStats : Base::SolverStats
    {
      IPMSolverStats()
      : Base::SolverStats()
      , delassus_decomposition_update_count(0)
      {
      }

      explicit IPMSolverStats(const int max_it)
      : it(0)
      , delassus_decomposition_update_count(0)
      {
        primal_feasibility.reserve(size_t(max_it));
        dual_feasibility.reserve(size_t(max_it));
        dual_feasibility_ncp.reserve(size_t(max_it));
        complementarity.reserve(size_t(max_it));
        mu.reserve(size_t(max_it));
      }

      void reset()
      {
        primal_feasibility.clear();
        dual_feasibility.clear();
        complementarity.clear();
        dual_feasibility_ncp.clear();
        mu.clear();
        it = 0;
        delassus_decomposition_update_count = 0;
      }

      size_t size() const
      {
        return primal_feasibility.size();
      }

      ///  \brief Number of total iterations.
      int it;

      ///  \brief Number of delassus decomposition updates.
      int delassus_decomposition_update_count;

      /// \brief History of primal feasibility values.
      std::vector<Scalar> primal_feasibility;

      /// \brief History of dual feasibility values.
      std::vector<Scalar> dual_feasibility;
      std::vector<Scalar> dual_feasibility_ncp;

      /// \brief History of complementarity values.
      std::vector<Scalar> complementarity;

      /// \brief History of mu (log barrier penalty) values.
      std::vector<Scalar> mu;
    };

    ///
    /// \brief Solve the constrained conic problem composed of problem data
    /// (G,g,cones) and starting from the initial guess.
    ///
    /// \param[in] G Symmetric PSD matrix representing the Delassus of the contact problem.
    /// \param[in] g Free contact acceleration or velocity associted with the contact problem.
    /// \param[in] constraint_models vector of constraint models
    /// \param[in] constraint_datas vector of constraint datas associated and up to date with
    /// constraint models
    /// \param[in] solve_ncp whether to solve the NCP (true) or CCP (false)
    /// \param[in] iterative_refinement_steps number of iterative refinement steps when solving PD
    /// Newton system
    /// \param[in] target_mu ?????
    ///
    /// \returns True if the problem has converged.
    template<
      typename DelassusDerived,
      typename VectorLikeConstraintDrift,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    bool solve(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLikeConstraintDrift> & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const boost::optional<RefConstVectorXs> primal_guess = boost::none,
      const boost::optional<RefConstVectorXs> dual_guess = boost::none,
      bool solve_ncp = true,
      bool stat_record = false,
      int iterative_refinement_steps = 0,
      Scalar target_mu = 1e-12,
      bool verbose = false);

    /// \returns the primal solution of the problem
    const VectorXs & getPrimalSolution() const
    {
      return x;
    }

    /// \returns the dual solution of the problem
    const VectorXs & getDualSolution() const
    {
      return z;
    }

    /// \returns the complementarity shift
    const VectorXs & getComplementarityShift() const
    {
      return saxce_corr;
    }

    /// \returns the stats of the solver (if solve was called with stat_record=true)
    const IPMSolverStats & getStats() const
    {
      return stats;
    }

  protected:
    using Base::problem_size;
    VectorXs v_constraint;
    VectorXs x, x_previous;
    VectorXs z, s, lambda;
    VectorXs delta_x, delta_z, delta_s;
    VectorXs delta_xi, delta_zi, delta_si;
    VectorXs ccp_primal_feas, ccp_primal_opt, ccp_compl_slackness, ncp_primal_opt;
    VectorXs tmp_vec_0, tmp_vec_1, tmp_vec_2, mehrotra_correction;
    VectorXs rhs_x, rhs_z, rhs_s;
    VectorXs rhs_x2, rhs_z2, rhs_s2;
    VectorXs saxce_corr;
    BarrierHessianTermVector barrier_hessian_terms;
    ScalingMatrixVector scaling_matrices;
    IPMSolverStats stats;
    bool v_constraint_evaluated = false;
    bool is_initialized = false;

    // algorithm options
    // TODO: public/setters
    /// \brief apply the corrector step at each iteration.
    bool corrector_step = true;

    /// \brief update the saxe correction at every iteration.
    // TODO: change to tol scheduling or numit for desaxce update
    bool always_update_saxe = false;

    /// \brief the barrier parameter / \mu
    Scalar barrier_parameter;

    bool second_order_saxe_correction = false;
    Scalar primal_opt_tol = 1e-6;

#ifdef PINOCCHIO_WITH_HPP_FCL
    using Base::timer;
#endif // PINOCCHIO_WITH_HPP_FCL

    template<typename DelassusDerived, typename VectorLikeConstraintDrift>
    const VectorXs & computeCurrentConstraintVelocity(
      DelassusOperatorBase<DelassusDerived> & delassus,
      const Eigen::MatrixBase<VectorLikeConstraintDrift> & g)
    {
      if (!v_constraint_evaluated)
      {
        delassus.applyOnTheRight(x, v_constraint);
        v_constraint += g;
        v_constraint_evaluated = true;
      }
      return v_constraint;
    }

    /// \brief Compute the primal feasibility vector (-Gx + s)
    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    void computePrimalFeasibilityVector(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

    /// \brief Compute the primal optimality vector (Px - G^T z + p + saxce)
    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    void computePrimalOptimalityVector(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

    /// \brief Compute the complementarity slackness vector (W^{-T} s ∘ W z)
    void computeComplementaritySlacknessVector();

    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator,
      typename VectorLikeInOut>
    static void normalizeConeVariables(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const Eigen::MatrixBase<VectorLikeInOut> & x);

    template<
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator,
      typename VectorLikeInOut>
    static void denormalizeConeVariables(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const Eigen::MatrixBase<VectorLikeInOut> & x);

    template<typename ConstraintModel, typename ConstraintModelAllocator>
    static void updateBarrierHessian(
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const ScalingMatrixVector & scaling_matrices,
      BarrierHessianTermVector & barrier_hessian_terms);

    /// \brief solve the primal-dual linear system
    //   Solve
    //
    //      [ 0     ]  +  [ G  T']   [ dx        ] = - [rhs_x]
    //      [ W'*ds ]  +  [ T  0 ]   [ W^{-1}*dz ] = - [rhs_z]
    //      lmbda o (dz + ds) = - [rhs_s]
    //      rhe right hand side is given by the object's member variables rhs_x,
    //      rhs_z, rhs_s. The solution is stored in the object's member variables
    //      dx, dz, ds. The matrix G is the delassus matrix, T is cone constrain
    //      matrix, W is the scaling matrix (member variable).
    //      The system is solved using a Schur complement approach to eliminate
    //      the dual variable z, after which a Cholesky factorization is used to solve for x.
    //      the result is dx, W*dz, W{-1}'*ds, which are stored in the object's
    //      delta_x, delta_s and delta_z member variables.
    //      Internally, the function uses the object's member variables tmp_vec_0,
    //      tmp_vec_1, tmp_vec_2, to store intermediate results without memory alloc.
    template<
      typename DelassusDerived,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    void solvePDSystem(
      const DelassusOperatorBase<DelassusDerived> & delassus,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      int iterative_refinement_steps);

    static void printIterationDetails(
      int iteration,
      Scalar primal_feas,
      Scalar primal_opt,
      Scalar gap,
      Scalar step_length,
      Scalar ncp_convergence);

    static void printIterationsHeader();

  }; // struct IPMConstraintSolverTpl

} // namespace pinocchio

#include "ipm-solver.hxx"

#endif // ifndef __pinocchio_algorithm_ip_solver_hpp__
