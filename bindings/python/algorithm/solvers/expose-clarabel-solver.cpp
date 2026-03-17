//
// Copyright (c) 2024-2025 INRIA
//

#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT

  #include "pinocchio/bindings/python/fwd.hpp"
  #include "pinocchio/algorithm/solvers/clarabel-solver.hpp"
  #include "pinocchio/algorithm/constraint-cholesky.hpp"
  #include "pinocchio/algorithm/delassus-operator-dense.hpp"
  #include "pinocchio/algorithm/delassus-operator-sparse.hpp"

  #include "pinocchio/bindings/python/utils/std-vector.hpp"
  #include "pinocchio/bindings/python/utils/macros.hpp"

  #include <boost/optional.hpp>

  #include <eigenpy/memory.hpp>
  #include <eigenpy/eigen-from-python.hpp>
  #include <eigenpy/eigen-to-python.hpp>

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef context::Scalar Scalar;
    static constexpr int Options = context::Options;
    typedef context::VectorXs VectorXs;
    typedef context::MatrixXs MatrixXs;

    typedef ClarabelConstraintSolverTpl<Scalar, Options> ClarabelConstraintSolver;
    typedef typename ClarabelConstraintSolver::ClarabelSolverSettings ClarabelSolverSettings;
    typedef typename ClarabelConstraintSolver::ClarabelSolverResult ClarabelSolverResult;

    typedef ConstraintSolverSettingsBaseTpl<Scalar> ConstraintSolverSettingsBase;
    typedef ConstraintSolverResultBaseTpl<Scalar> ConstraintSolverResultBase;
    typedef ConstraintSolverStatsBaseTpl<Scalar> ConstraintSolverStatsBase;
    typedef ConstraintSolverBaseTpl<Scalar> ConstraintSolverBase;

    typedef ContactCholeskyDecompositionTpl<Scalar, Options> ContactCholeskyDecomposition;

    // ============================================================================
    // Expose ClarabelSolverSettings (inheriting from base)
    // ============================================================================

    void exposeClarabelSolverSettings()
    {
      bp::class_<ClarabelSolverSettings, bp::bases<ConstraintSolverSettingsBase>>(
        "ClarabelSolverSettings", "Settings for the clarabel constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor with default settings."))

        // Clarabel specific settings (base class properties are inherited)
        .def_readwrite(
          "max_ncp_loops", &ClarabelSolverSettings::max_ncp_loops,
          "Maximum number of calls to the clarabel solver when "
          "solving the NCP problem.")
        .def_readwrite(
          "tol_ktratio", &ClarabelSolverSettings::tol_ktratio,
          "Tolerance for KKT conditions (Clarabel-specific).")
        .def_readwrite(
          "verbose", &ClarabelSolverSettings::verbose,
          "Enable verbose output (Clarabel-specific).");
    }

    // ============================================================================
    // Expose ClarabelSolverResult (inheriting from base)
    // ============================================================================

    // Wrapper functions for retrieve methods
    static void retrievePrimalSolution_wrapper(
      const ClarabelSolverResult & solution, Eigen::Ref<VectorXs> primal_solution)
    {
      solution.retrievePrimalSolution(primal_solution);
    }

    static void retrieveDualSolution_wrapper(
      const ClarabelSolverResult & solution, Eigen::Ref<VectorXs> dual_solution)
    {
      solution.retrieveDualSolution(dual_solution);
    }

    static void retrieveDesaxceTerm_wrapper(
      const ClarabelSolverResult & solution, Eigen::Ref<VectorXs> desaxce_term)
    {
      solution.retrieveDesaxceTerm(desaxce_term);
    }

    void exposeClarabelSolverResult()
    {
      bp::class_<ClarabelSolverResult, bp::bases<ConstraintSolverResultBase>>(
        "ClarabelSolverResult", "Solution of the Clarabel constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor."))

        .PINOCCHIO_ADD_PROPERTY_READONLY(ClarabelSolverResult, problem_size, "Problem size")

        .def(
          "reset",
          static_cast<void (ClarabelSolverResult::*)(std::size_t)>(&ClarabelSolverResult::reset),
          (bp::arg("self"), bp::arg("problem_size") = 0), "Reset the result")
        .def(
          "resize", &ClarabelSolverResult::resize, bp::args("self", "problem_size"),
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
    // Expose ClarabelConstraintSolver (inheriting from base)
    // ============================================================================

  #ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
    template<
      typename DelassusDerived,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    static bool solve_clarabel_wrapper(
      ClarabelConstraintSolver & solver,
      DelassusDerived & delassus,
      const VectorXs & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const ClarabelSolverSettings & settings,
      ClarabelSolverResult & result)
    {
      return solver.solve(delassus, g, constraint_models, constraint_datas, settings, result);
    }
  #endif

    template<typename Solver>
    struct ClarabelSolveMethodExposer
    {
      ClarabelSolveMethodExposer(bp::class_<Solver, bp::bases<ConstraintSolverBase>> & class_)
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
        typedef std::allocator<ConstraintModel> ConstraintModelAllocator;
        typedef std::allocator<ConstraintData> ConstraintDataAllocator;

        PINOCCHIO_UNUSED_VARIABLE(ptr);

  #ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
        class_
          .def(
            "solve",
            solve_clarabel_wrapper<
              ContactCholeskyDecomposition::DelassusCholeskyExpression, ConstraintModel,
              ConstraintModelAllocator, ConstraintData, ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.")
          .def(
            "solve",
            solve_clarabel_wrapper<
              context::DelassusOperatorDense, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.")
          .def(
            "solve",
            solve_clarabel_wrapper<
              context::DelassusOperatorSparse, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.");

    #ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
        {
          typedef Eigen::AccelerateLLT<context::SparseMatrix> AccelerateLLT;
          typedef DelassusOperatorSparseTpl<context::Scalar, Options, AccelerateLLT>
            DelassusOperatorSparseAccelerate;
          class_.def(
            "solve",
            solve_clarabel_wrapper<
              DelassusOperatorSparseAccelerate, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.");
        }
    #endif // ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
  #endif   // ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
      }

      void run(boost::blank * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

      bp::class_<Solver, bp::bases<ConstraintSolverBase>> & class_;
    };

    void exposeClarabelConstraintSolver()
    {
  #ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE

      // Expose Settings, Solution, Stats (they inherit from base)
      exposeClarabelSolverSettings();
      exposeClarabelSolverResult();

      // Expose the solver itself (inherits from base)
      bp::class_<ClarabelConstraintSolver, bp::bases<ConstraintSolverBase>> cl(
        "ClarabelConstraintSolver",
        "Alternating Direction Method of Multipliers (Clarabel) solver for contact dynamics.",
        bp::init<std::size_t>(
          bp::args("self", "problem_size"), "Constructor with problem dimension."));

      cl.def(
          "isValid", &ClarabelConstraintSolver::isValid, bp::arg("self"),
          "Check if the solver is in a valid state (has solved a constraint problem)")
        .def(
          "reset", &ClarabelConstraintSolver::reset, bp::arg("self"),
          "Reset the solver to initial state");

      // Expose solve methods for different constraint models
      ClarabelSolveMethodExposer<ClarabelConstraintSolver> solve_exposer(cl);
      solve_exposer.run(static_cast<context::ConstraintModel *>(nullptr));

  #endif
    }

  } // namespace python
} // namespace pinocchio

#endif // PINOCCHIO_WITH_CLARABEL_SUPPORT
