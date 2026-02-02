//
// Copyright (c) 2018 CNRS
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/multibody/liegroup.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<typename LieGroupCollection>
  struct LieGroupGenericTpl;

} // namespace pinocchio
