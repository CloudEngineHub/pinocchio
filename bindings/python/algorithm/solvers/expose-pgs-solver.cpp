//
// Copyright (c) 2022-2025 INRIA
//

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/algorithm/solvers/pgs-solver.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"

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

    typedef PGSConstraintSolverTpl<Scalar> PGSSolver;
    typedef typename PGSSolver::PGSSolverSettings PGSSolverSettings;
    typedef typename PGSSolver::PGSSolverResult PGSSolverResult;
    typedef typename PGSSolver::PGSSolverStats PGSSolverStats;

    typedef ConstraintSolverSettingsBaseTpl<Scalar> ConstraintSolverSettingsBase;
    typedef ConstraintSolverResultBaseTpl<Scalar> ConstraintSolverResultBase;
    typedef ConstraintSolverStatsBaseTpl<Scalar> ConstraintSolverStatsBase;
    typedef ConstraintSolverBaseTpl<Scalar> ConstraintSolverBase;

    // ============================================================================
    // Expose PGSSolverSettings (inheriting from base)
    // ============================================================================

    void exposePGSSolverSettings()
    {
      bp::class_<PGSSolverSettings, bp::bases<ConstraintSolverSettingsBase>>(
        "PGSSolverSettings", "Settings for the PGS constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor with default settings."))

        // PGS specific settings (base class properties are inherited)
        .PINOCCHIO_ADD_PROPERTY(
          PGSSolverSettings, over_relaxation,
          "Over-relaxation parameter (should be in ]0,2[, default 1)");
    }

    // ============================================================================
    // Expose PGSSolverResult (inheriting from base)
    // ============================================================================

    // Wrapper functions for retrieve methods
    static void
    retrievePrimalSolution_pgs_wrapper(const PGSSolverResult & solution, VectorXs & primal_solution)
    {
      solution.retrievePrimalSolution(primal_solution);
    }

    static void
    retrieveDualSolution_pgs_wrapper(const PGSSolverResult & solution, VectorXs & dual_solution)
    {
      solution.retrieveDualSolution(dual_solution);
    }

    void exposePGSSolverResult()
    {
      bp::class_<PGSSolverResult, bp::bases<ConstraintSolverResultBase>>(
        "PGSSolverResult", "Solution of the PGS constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor."))

        // PGS specific properties (base class properties are inherited)
        .PINOCCHIO_ADD_PROPERTY_READONLY(PGSSolverResult, problem_size, "Problem size")

        .def(
          "reset", static_cast<void (PGSSolverResult::*)(std::size_t)>(&PGSSolverResult::reset),
          (bp::arg("self"), bp::arg("problem_size") = 0), "Reset the result")
        .def(
          "resize", &PGSSolverResult::resize, bp::args("self", "problem_size"),
          "Resize solution vectors")

        // Retrieve methods
        .def(
          "retrievePrimalSolution", retrievePrimalSolution_pgs_wrapper,
          bp::args("self", "primal_solution"),
          "Retrieve the primal solution (copy x into primal_solution)")
        .def(
          "retrieveDualSolution", retrieveDualSolution_pgs_wrapper,
          bp::args("self", "dual_solution"),
          "Retrieve the dual solution (copy y into dual_solution)");
    }

    // ============================================================================
    // Expose PGSSolverStats (inheriting from base)
    // ============================================================================

    void exposePGSSolverStats()
    {
      bp::class_<PGSSolverStats, bp::bases<ConstraintSolverStatsBase>>(
        "PGSSolverStats", "Per-iteration statistics of the PGS constraint solver.",
        bp::init<>(bp::arg("self"), "Default constructor."))
        .def(
          bp::init<std::size_t>(
            bp::args("self", "max_iterations"), "Constructor with maximum iterations."));
      // Note: No PGS-specific stats beyond base class
    }

    // ============================================================================
    // Expose PGSConstraintSolver (inheriting from base)
    // ============================================================================

#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
    template<
      typename DelassusDerived,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    static bool solve_pgs_wrapper(
      PGSSolver & solver,
      const DelassusDerived & delassus,
      const VectorXs & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const PGSSolverSettings & settings,
      PGSSolverResult & result)
    {
      return solver.solve(delassus, g, constraint_models, constraint_datas, settings, result);
    }
#endif

    template<typename Solver>
    struct PGSSolveMethodExposer
    {
      PGSSolveMethodExposer(bp::class_<Solver, bp::bases<ConstraintSolverBase>> & class_)
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

#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
        class_
          .def(
            "solve",
            solve_pgs_wrapper<
              context::MatrixXs, ConstraintModel, ConstraintModelAllocator, ConstraintData,
              ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.")
          .def(
            "solve",
            solve_pgs_wrapper<
              context::DelassusOperatorDense, ConstraintModel, ConstraintModelAllocator,
              ConstraintData, ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.")
          .def(
            "solve",
            solve_pgs_wrapper<
              ContactCholeskyDecomposition::DelassusCholeskyExpression, ConstraintModel,
              ConstraintModelAllocator, ConstraintData, ConstraintDataAllocator>,
            bp::args(
              "self", "delassus", "g", "constraint_models", "constraint_datas", "settings",
              "result"),
            "Solve the constrained conic problem with given settings and result.");
#endif // ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
      }

      void run(boost::blank * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

      bp::class_<Solver, bp::bases<ConstraintSolverBase>> & class_;
    };

    void exposePGSConstraintSolver()
    {
#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE

      // Expose Settings, Solution, Stats (they inherit from base)
      exposePGSSolverSettings();
      exposePGSSolverResult();
      exposePGSSolverStats();

      // Expose the solver itself (inherits from base)
      bp::class_<PGSSolver, bp::bases<ConstraintSolverBase>> cl(
        "PGSConstraintSolver", "Projected Gauss-Seidel (PGS) solver for contact dynamics.",
        bp::init<std::size_t>(
          bp::args("self", "problem_size"), "Constructor with problem dimension."));

      cl.PINOCCHIO_ADD_PROPERTY_READONLY(PGSSolver, stats, "Access the statistics of the solver")
        .def(
          "isValid", &PGSSolver::isValid, bp::arg("self"),
          "Check if the solver is in a valid state (has solved a constraint problem)")
        .def("reset", &PGSSolver::reset, bp::arg("self"), "Reset the solver to initial state");

      // Expose solve methods for different constraint models
      PGSSolveMethodExposer<PGSSolver> solve_exposer(cl);
      solve_exposer.run(static_cast<context::ConstraintModel *>(nullptr));

#endif
    }

  } // namespace python
} // namespace pinocchio
