
#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/multibody/liegroup.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<
    typename Scalar,
    int Options = context::Options,
    template<typename, int> class LieGroupCollectionTpl = LieGroupCollectionDefaultTpl>
  struct CartesianProductOperationVariantTpl;

  typedef CartesianProductOperationVariantTpl<
    context::Scalar,
    context::Options,
    LieGroupCollectionDefaultTpl>
    CartesianProductOperationVariant;
} // namespace pinocchio
