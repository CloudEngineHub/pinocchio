//
// Copyright (c) 2022-2024 INRIA
//

#include "pinocchio/bindings/python/algorithm/algorithms.hpp"

namespace pinocchio
{
  namespace python
  {

    // Forward declarations
    void exposePGSConstraintSolver();
    void exposeADMMConstraintSolver();
#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT
    void exposeClarabelConstraintSolver();
#endif

    void exposeConstraintSolvers()
    {
      exposePGSConstraintSolver();
      exposeADMMConstraintSolver();
#ifdef PINOCCHIO_WITH_CLARABEL_SUPPORT
      exposeClarabelConstraintSolver();
#endif
    }

  } // namespace python
} // namespace pinocchio
