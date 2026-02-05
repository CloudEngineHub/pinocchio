//
// Copyright (c) 2025 INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<typename T, typename Enable = void>
  struct sizeInBytesImpl;
} // namespace pinocchio
