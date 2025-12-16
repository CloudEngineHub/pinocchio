//
// Copyright (c) 2024-2025 INRIA
//

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/algorithm/solvers/constraint-solver-base.hpp"

#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef context::Scalar Scalar;

    typedef ConstraintSolverSettingsBaseTpl<Scalar> ConstraintSolverSettingsBase;
    typedef ConstraintSolverSolutionBaseTpl<Scalar> ConstraintSolverSolutionBase;
    typedef ConstraintSolverStatsBaseTpl<Scalar> ConstraintSolverStatsBase;
    typedef ConstraintSolverBaseTpl<Scalar> ConstraintSolverBase;

    // ============================================================================
    // Expose ConstraintSolverSettingsBase
    // ============================================================================

    void exposeConstraintSolverSettingsBase()
    {
      bp::class_<ConstraintSolverSettingsBase>(
        "ConstraintSolverSettingsBase", "Base class for constraint solver settings.",
        bp::no_init) // Abstract base class, no direct construction

        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, max_iterations, "Maximum number of iterations")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, absolute_tol_feasibility,
          "Absolute tolerance on primal/dual feasibility")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, relative_tol_feasibility,
          "Relative tolerance on primal/dual feasibility")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, absolute_tol_complementarity,
          "Absolute tolerance on complementarity")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, relative_tol_complementarity,
          "Relative tolerance on complementarity")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, solve_ncp, "Whether to solve NCP (true) or CCP (false)")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, measure_timings, "Whether to measure timings")
        .PINOCCHIO_ADD_PROPERTY(
          ConstraintSolverSettingsBase, stat_record, "Whether to record per-iteration statistics")

        .def(
          "checkValidity", &ConstraintSolverSettingsBase::checkValidity, bp::arg("self"),
          "Check if settings are valid (throws if invalid)");
    }

    // ============================================================================
    // Expose ConstraintSolverSolutionBase
    // ============================================================================

    void exposeConstraintSolverSolutionBase()
    {
      bp::class_<ConstraintSolverSolutionBase>(
        "ConstraintSolverSolutionBase", "Base class for constraint solver solution.",
        bp::no_init) // Abstract base class, no direct construction

        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverSolutionBase, iterations, "Number of iterations")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverSolutionBase, converged, "Whether solver converged")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverSolutionBase, primal_feasibility, "Final primal feasibility")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverSolutionBase, dual_feasibility, "Final dual feasibility")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverSolutionBase, complementarity, "Final complementarity")

        .def("reset", &ConstraintSolverSolutionBase::reset, bp::arg("self"), "Reset the solution");
    }

    // ============================================================================
    // Expose ConstraintSolverStatsBase
    // ============================================================================

    void exposeConstraintSolverStatsBase()
    {
      bp::class_<ConstraintSolverStatsBase>(
        "ConstraintSolverStatsBase", "Base class for constraint solver per-iteration statistics.",
        bp::no_init) // Abstract base class, no direct construction

        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverStatsBase, iterations, "Number of iterations")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverStatsBase, primal_feasibility, "History of primal feasibility")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverStatsBase, dual_feasibility, "History of dual feasibility")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverStatsBase, dual_feasibility_ncp, "History of NCP dual feasibility")
        .PINOCCHIO_ADD_PROPERTY_READONLY(
          ConstraintSolverStatsBase, complementarity, "History of complementarity")

        .def("reset", &ConstraintSolverStatsBase::reset, bp::arg("self"), "Reset the statistics")
        .def(
          "reserve", &ConstraintSolverStatsBase::reserve, bp::args("self", "max_iterations"),
          "Reserve storage for max_iterations")
        .def(
          "size", &ConstraintSolverStatsBase::size, bp::arg("self"),
          "Size of the statistics (number of recorded iterations)");
    }

    // ============================================================================
    // Expose ConstraintSolverBase
    // ============================================================================

    void exposeConstraintSolverBase()
    {
      bp::class_<ConstraintSolverBase>(
        "ConstraintSolverBase", "Base class for constraint solvers.",
        bp::no_init) // Abstract base class, no direct construction

#ifdef PINOCCHIO_WITH_HPP_FCL
        .def(
          "getCPUTimes", &ConstraintSolverBase::getCPUTimes, bp::arg("self"),
          "Get CPU times for the solve")
#endif
        ;
    }

    // ============================================================================
    // Main exposure function
    // ============================================================================

    void exposeConstraintSolverBases()
    {
#ifdef PINOCCHIO_PYTHON_PLAIN_SCALAR_TYPE
      exposeConstraintSolverBase();
      exposeConstraintSolverSettingsBase();
      exposeConstraintSolverSolutionBase();
      exposeConstraintSolverStatsBase();
#endif
    }

  } // namespace python
} // namespace pinocchio
