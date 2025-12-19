//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_std_array_hpp__
#define __pinocchio_utils_std_array_hpp__

#include <array>

#include "pinocchio/utils/size-in-bytes.hpp"

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

#endif // __pinocchio_utils_std_array_hpp__
