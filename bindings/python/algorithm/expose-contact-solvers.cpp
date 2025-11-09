//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"

namespace pinocchio
{
  namespace python
  {

    // Forward declaration
    void exposePGSContactSolver();
    void exposeADMMContactSolver();
    void exposeIPMConstraintSolver();

    void exposeContactSolvers()
    {
      exposePGSContactSolver();
      exposeADMMContactSolver();
      exposeIPMConstraintSolver();
    }

  } // namespace python
} // namespace pinocchio
