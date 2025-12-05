//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"

namespace pinocchio
{
  namespace python
  {

    // Forward declaration
    void exposePGSConstraintSolver();
    void exposeADMMConstraintSolver();
    void exposeIPMConstraintSolver();

    void exposeContactSolvers()
    {
      exposePGSConstraintSolver();
      exposeADMMConstraintSolver();
      exposeIPMConstraintSolver();
    }

  } // namespace python
} // namespace pinocchio
