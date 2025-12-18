//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_size_in_bytes_hpp__
#define __pinocchio_utils_size_in_bytes_hpp__

#include <cstddef>
#include <type_traits>
#include <utility> // for std::declval

namespace pinocchio
{
  namespace helper
  {
    template<typename T, typename = void>
    struct has_method_sizeInBytes : std::false_type
    {
    };

    template<typename T>
    struct has_method_sizeInBytes<T, std::void_t<decltype(std::declval<T &>().sizeInBytes())>>
    : std::true_type
    {
    };

    template<typename T>
    inline constexpr bool has_method_sizeInBytes_v = has_method_sizeInBytes<T>::value;
  } // namespace helper

  template<typename T>
  struct sizeInBytesImpl
  {
    static std::size_t run(const T & value);
  };

  template<typename T>
  std::size_t sizeInBytes(const T & value)
  {
    if constexpr (helper::has_method_sizeInBytes_v<T>)
    {
      return value.sizeInBytes();
    }
    else
      return sizeInBytesImpl<T>::run(value);
  }

} // namespace pinocchio

#endif // __pinocchio_utils_size_in_bytes_hpp__
