//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"

namespace pinocchio
{
  namespace python
  {

    // Forward declarations
    void exposeConstraintSolverBases(); // Expose base classes first
    void exposePGSConstraintSolver();
    void exposeADMMConstraintSolver();

    void exposeConstraintSolvers()
    {
      exposeConstraintSolverBases(); // Must be called first!
      exposePGSConstraintSolver();
      exposeADMMConstraintSolver();
    }

  } // namespace python
} // namespace pinocchio
