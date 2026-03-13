//
// Copyright (c) 2024 INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include <memory>
#endif // PINOCCHIO_LSP

namespace coal
{
  class MeshLoader;
  typedef std::shared_ptr<MeshLoader> MeshLoaderPtr;
} // namespace coal
