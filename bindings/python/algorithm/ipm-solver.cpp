//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/bindings/python/fwd.hpp"

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/algorithm/ipm-solver.hpp"
#include "pinocchio/algorithm/delassus-operator-dense.hpp"

#include "pinocchio/bindings/python/algorithm/contact-solver-base.hpp"
#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"
#include <eigenpy/eigen-from-python.hpp>

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    typedef context::Scalar Scalar;
    typedef context::VectorXs VectorXs;
    typedef IPMConstraintSolverTpl<Scalar> Solver;
    typedef Solver::IPMSolverStats SolverStats;
    typedef typename Solver::RefConstVectorXs RefConstVectorXs;

    template<
      typename DelassusDerived,
      typename ConstraintModel,
      typename ConstraintModelAllocator,
      typename ConstraintData,
      typename ConstraintDataAllocator>
    static bool solve_wrapper(
      Solver & solver,
      DelassusDerived & delassus,
      const VectorXs & g,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const boost::optional<RefConstVectorXs> primal_guess = boost::none,
      const boost::optional<RefConstVectorXs> dual_guess = boost::none,
      const bool solve_ncp = true,
      const bool stat_record = false,
      int iterative_refinenment_steps = 0,
      Scalar target_mu = 1e-12,
      bool verbose = false)
    {
      return solver.solve(
        delassus, g, constraint_models, constraint_datas, primal_guess, dual_guess, //
        solve_ncp, stat_record, iterative_refinenment_steps, target_mu, verbose);
    }

    template<typename Solver>
    struct SolveMethodExposer
    {
      SolveMethodExposer(bp::class_<Solver> & class_)
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
        PINOCCHIO_UNUSED_VARIABLE(ptr);
        typedef typename traits<ConstraintModel>::ConstraintData ConstraintData;
        typedef Eigen::aligned_allocator<ConstraintModel> ConstraintModelAllocator;
        typedef Eigen::aligned_allocator<ConstraintData> ConstraintDataAllocator;

        class_.def(
          "solve",
          solve_wrapper<
            context::DelassusOperatorDense, ConstraintModel, ConstraintModelAllocator,
            ConstraintData, ConstraintDataAllocator>,
          (bp::args("self", "delassus", "g", "constraint_models", "constraint_datas"),
           bp::arg("primal_guess") = boost::none, bp::arg("dual_guess") = boost::none,
           bp::arg("solve_ncp") = true, bp::arg("stat_record") = false,
           bp::arg("iterative_refinement_steps") = 0, bp::arg("target_mu") = 1e-12,
           bp::arg("verbose") = false),
          "Solve the constrained conic problem composed of problem data (G,g,cones) and starting "
          "from the initial guess.");
      }

      void run(boost::blank * ptr = 0)
      {
        PINOCCHIO_UNUSED_VARIABLE(ptr);
      }

      bp::class_<Solver> & class_;
    };

    template<typename ConstraintModel>
    static void expose_solve(bp::class_<Solver> & class_)
    {
      SolveMethodExposer<Solver> expose(class_);
      expose.run(static_cast<ConstraintModel *>(nullptr));
    }

    void exposeIPMConstraintSolver()
    {
      // Expose constructor
      bp::class_<Solver> class_(
        "IPMConstraintSolver", "Interior Point Method solver for constraint dynamics.",
        bp::init<int>(bp::args("self", "problem_dim"), "Default constructor."));

      // Expose helpers
      class_.def(ConstraintSolverBasePythonVisitor<Solver>())
        .def(
          "getPrimalSolution", &Solver::getPrimalSolution, bp::arg("self"),
          "Returns the primal solution of the problem.", bp::return_internal_reference<>())
        .def(
          "getDualSolution", &Solver::getDualSolution, bp::arg("self"),
          "Returns the dual solution of the problem.", bp::return_internal_reference<>())
        .def("getStats", &Solver::getStats, bp::arg("self"), bp::return_internal_reference<>());

      // Expose solve
      {
        expose_solve<context::ConstraintModel>(class_);
      }

      // Expose stats
      {
        bp::class_<SolverStats>(
          "SolverStats", "",
          bp::init<int>((bp::arg("self"), bp::arg("max_it")), "Default constructor"))
          .def("reset", &SolverStats::reset, bp::arg("self"), "Reset the stasts.")
          .def(
            "size", &SolverStats::size, bp::arg("self"),
            "Size of the vectors stored in the structure.")

          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, primal_feasibility, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, dual_feasibility, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, dual_feasibility_ncp, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, complementarity, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, mu, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(SolverStats, delassus_decomposition_update_count, "")
          .PINOCCHIO_ADD_PROPERTY_READONLY(
            SolverStats, it, "Number of iterations performed by the algorithm.");
      }
    }

  } // namespace python
} // namespace pinocchio
