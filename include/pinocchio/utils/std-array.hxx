//
// Copyright (c) 2025 INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/utils/size-in-bytes.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  template<typename T, std::size_t N>
  struct sizeInBytesImpl<std::array<T, N>>
  {
    static std::size_t run(const std::array<T, N> & array)
    {
      std::size_t size_value = 0;
      for (const auto & elt : array)
      {
        size_value += sizeInBytes(elt);
      }
      return size_value;
    }
  }; // sizeInBytesImpl

} // namespace pinocchio
