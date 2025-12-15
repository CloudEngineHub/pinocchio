//
// Copyright (c) 2024-2025 INRIA
//

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/algorithm/solvers/admm-solver.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include "pinocchio/algorithm/delassus-operator-sparse.hpp"

#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

#include <eigenpy/memory.hpp>
#include <eigenpy/eigen-from-python.hpp>
#include <eigenpy/eigen-to-python.hpp>

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef context::Scalar Scalar;
    typedef context::VectorXs VectorXs;
    typedef context::MatrixXs MatrixXs;

    typedef ADMMConstraintSolverTpl<Scalar> ADMMSolver;
    typedef typename ADMMSolver::ADMMSolverSettings ADMMSolverSettings;
    typedef typename ADMMSolver::ADMMSolverSolution ADMMSolverSolution;
    typedef typename ADMMSolver::ADMMSolverStats ADMMSolverStats;

    typedef ConstraintSolverSettingsBaseTpl<Scalar> ConstraintSolverSettingsBase;
    typedef ConstraintSolverSolutionBaseTpl<Scalar> ConstraintSolverSolutionBase;
    typedef ConstraintSolverStatsBaseTpl<Scalar> ConstraintSolverStatsBase;
    typedef ConstraintSolverBaseTpl<Scalar> ConstraintSolverBase;

    typedef ContactCholeskyDecompositionTpl<Scalar, context::Options> ContactCholeskyDecomposition;

    // ============================================================================
    // Expose ADMM Enums
    // ============================================================================

    void exposeADMMEnums()
    {
      if (!eigenpy::register_symbolic_link_to_registered_type<::pinocchio::ADMMUpdateRule>())
      {
        bp::enum_<::pinocchio::ADMMUpdateRule>("ADMMUpdateRule")
          .value("SPECTRAL", ::pinocchio::ADMMUpdateRule::SPECTRAL)
          .value("OSQP", ::pinocchio::ADMMUpdateRule::OSQP)
          .value("LINEAR", ::pinocchio::ADMMUpdateRule::LINEAR)
          .value("CONSTANT", ::pinocchio::ADMMUpdateRule::CONSTANT)
          .export_values();
      }

      if (!eigenpy::register_symbolic_link_to_registered_type<::pinocchio::ADMMProximalRule>())
      {
        bp::enum_<::pinocchio::ADMMProximalRule>("ADMMProximalRule")
          .value("MANUAL", ::pinocchio::ADMMProximalRule::MANUAL)
          .value("AUTOMATIC", ::pinocchio::ADMMProximalRule::AUTOMATIC)
          .export_values();
      }
    }

    // ============================================================================
    // Expose ADMMSolverSettings (inheriting from base)
    // ============================================================================

    void exposeADMMSolverSettings()
    {
      bp::class_<ADMMSolverSettings, bp::bases<ConstraintSolverSettingsBase>>(
        "ADMMSolverSettings", "Settings for the ADMM constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor with default settings."))

        // ADMM specific settings (base class properties are inherited)
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, admm_update_rule, "ADMM update rule")
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, admm_proximal_rule, "ADMM proximal rule")
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, mu_prox, "Proximal penalty parameter")
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, tau_prox, "Proximal scaling factor")
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, tau, "ADMM penalty scaling factor")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, ratio_primal_dual, "Primal/dual ratio threshold for rho update")
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, dual_momentum, "Dual variable momentum")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, rho_update_ratio, "Ratio threshold for rho update")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, rho_min_update_frequency, "Minimum frequency for rho updates")
        .PINOCCHIO_ADD_PROPERTY(ADMMSolverSettings, rho_momentum, "Momentum on rho updates")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, spectral_rho_power_init, "Initial rho power for SPECTRAL rule")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, spectral_rho_power_factor, "Rho power factor for SPECTRAL rule")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, linear_update_rule_factor, "Factor for LINEAR update rule")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, lanczos_size,
          "Size of Lanczos decomposition for eigenvalue estimation")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, max_delassus_decomposition_updates,
          "Maximum number of Delassus decomposition updates")
        .PINOCCHIO_ADD_PROPERTY(
          ADMMSolverSettings, anderson_capacity, "Anderson acceleration capacity");
    }

    // ============================================================================
    // Expose ADMMSolverSolution (inheriting from base)
    // ============================================================================

    // Wrapper functions for retrieve methods
    static void
    retrievePrimalSolution_wrapper(const ADMMSolverSolution & solution, VectorXs & primal_solution)
    {
      solution.retrievePrimalSolution(primal_solution);
    }

    static void
    retrieveDualSolution_wrapper(const ADMMSolverSolution & solution, VectorXs & dual_solution)
    {
      solution.retrieveDualSolution(dual_solution);
    }

    static void
    retrieveDesaxceTerm_wrapper(const ADMMSolverSolution & solution, VectorXs & desaxce_term)
    {
      solution.retrieveDesaxceTerm(desaxce_term);
    }

    void exposeADMMSolverSolution()
    {
      bp::class_<ADMMSolverSolution, bp::bases<ConstraintSolverSolutionBase>>(
        "ADMMSolverSolution", "Solution of the ADMM constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor."))

        // ADMM specific properties (base class properties are inherited)
        .PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolverSolution, problem_size, "Problem size")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ADMMSolverSolution, delassus_decomposition_update_count,
          "Number of Delassus decomposition updates")
        .PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolverSolution, rho, "Final rho value")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ADMMSolverSolution, spectral_rho_power, "Final spectral rho power")
        .PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolverSolution, mu_prox, "Final proximal parameter")

        .def(
          "resize", &ADMMSolverSolution::resize, bp::args("self", "problem_size"),
          "Resize solution vectors")

        // Retrieve methods
        .def(
          "retrievePrimalSolution", retrievePrimalSolution_wrapper,
          bp::args("self", "primal_solution"),
          "Retrieve the primal solution (copy y into primal_solution)")
        .def(
          "retrieveDualSolution", retrieveDualSolution_wrapper, bp::args("self", "dual_solution"),
          "Retrieve the dual solution (copy z into dual_solution)")
        .def(
          "retrieveDesaxceTerm", retrieveDesaxceTerm_wrapper, bp::args("self", "desaxce_term"),
          "Retrieve the DeSaxce correction term (copy desaxce into desaxce_term)");
    }

    // ============================================================================
    // Expose ADMMSolverStats (inheriting from base)
    // ============================================================================

    void exposeADMMSolverStats()
    {
      bp::class_<ADMMSolverStats, bp::bases<ConstraintSolverStatsBase>>(
        "ADMMSolverStats", "Per-iteration statistics of the ADMM constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor."))
        .def(bp::init<std::size_t>(
          bp::args("self", "max_iterations"), "Constructor with maximum iterations."))

        // ADMM specific properties (base class properties are inherited)
        .PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolverStats, rho, "History of rho values")
        .PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolverStats, mu_prox, "History of mu_prox values")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ADMMSolverStats, anderson_size, "History of Anderson acceleration size")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ADMMSolverStats, linear_system_residual, "History of linear system residuals")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ADMMSolverStats, linear_system_consistency, "History of linear system consistency")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ADMMSolverStats, delassus_decomposition_update_count,
          "Number of Delassus decomposition updates");
    }

    // ============================================================================
    // Expose ADMMConstraintSolver (inheriting from base)
    // ============================================================================

#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
    template<
      typename DelassusDerived,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    static bool solve_admm_wrapper(
      ADMMSolver & solver,
      DelassusDerived & delassus,
      const VectorXs & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const ADMMSolverSettings & settings)
    {
      return solver.solve(delassus, g, constraint_models, constraint_datas, settings);
    }
#endif

    template<typename Solver>
    struct ADMMSolveMethodExposer
    {
      ADMMSolveMethodExposer(bp::class_<Solver, bp::bases<ConstraintSolverBase>> & class_)
      : class_(class_)
      {
      }

      template<class T>
      void operator()(T)
      {
        run(static_cast<typename T::type *>(nullptr));
      }

      template<typename ConstraintModel>
      void run(ConstraintModelBase<ConstraintModel> * ptr = 0)
      {
        using ConstraintData = typename traits<ConstraintModel>::ConstraintData;
        typedef Eigen::aligned_allocator<ConstraintModel> ConstraintModelAllocator;
        typedef Eigen::aligned_allocator<ConstraintData> ConstraintDataAllocator;

        PINOCCHIO_UNUSED_VARIABLE(ptr);

        class_
          .def(
            "solve",
            solve_admm_wrapper<
              ContactCholeskyDecomposition::DelassusCholeskyExpression, ConstraintModel,
              ConstraintModelAllocator, ConstraintData, ConstraintDataAllocator>,
            bp::args("self", "delassus", "g", "constraint_models", "constraint_datas", "settings"),
            "Solve the constrained conic problem with given settings.")
          .def(
            "solve",
            solve_admm_wrapper<
              context::DelassusOperatorDense, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args("self", "delassus", "g", "constraint_models", "constraint_datas", "settings"),
            "Solve the constrained conic problem with given settings.")
          .def(
            "solve",
            solve_admm_wrapper<
              context::DelassusOperatorSparse, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args("self", "delassus", "g", "constraint_models", "constraint_datas", "settings"),
            "Solve the constrained conic problem with given settings.");

#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
        {
          typedef Eigen::AccelerateLLT<context::SparseMatrix> AccelerateLLT;
          typedef DelassusOperatorSparseTpl<context::Scalar, context::Options, AccelerateLLT>
            DelassusOperatorSparseAccelerate;
          class_.def(
            "solve",
            solve_admm_wrapper<
              DelassusOperatorSparseAccelerate, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args("self", "delassus", "g", "constraint_models", "constraint_datas", "settings"),
            "Solve the constrained conic problem with given settings.");
        }
#endif
      }

      void run(boost::blank * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

      bp::class_<Solver, bp::bases<ConstraintSolverBase>> & class_;
    };

    void exposeADMMConstraintSolver()
    {
#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE

      // Expose enums first
      exposeADMMEnums();

      // Expose Settings, Solution, Stats (they inherit from base)
      exposeADMMSolverSettings();
      exposeADMMSolverSolution();
      exposeADMMSolverStats();

      // Expose the solver itself (inherits from base)
      bp::class_<ADMMSolver, bp::bases<ConstraintSolverBase>> cl(
        "ADMMConstraintSolver",
        "Alternating Direction Method of Multipliers (ADMM) solver for contact dynamics.",
        bp::init<std::size_t>(
          bp::args("self", "problem_size"), "Constructor with problem dimension."));

      cl.PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolver, solution, "Access the solution of the solver")
        .PINOCCHIO_ADD_PROPERTY_READONLY(ADMMSolver, stats, "Access the statistics of the solver");

      // Expose solve methods for different constraint models
      ADMMSolveMethodExposer<ADMMSolver> solve_exposer(cl);
      solve_exposer.run(static_cast<context::ConstraintModel *>(nullptr));

#endif
    }

  } // namespace python
} // namespace pinocchio
