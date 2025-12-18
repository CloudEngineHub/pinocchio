//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_size_in_bytes_hpp__
#define __pinocchio_utils_size_in_bytes_hpp__

#include <type_traits>
#include <utility> // for std::declval

namespace pinocchio
{
  namespace helper
  {
    template<typename template<>, typename = void>
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

} // namespace pinocchio

#endif // __pinocchio_utils_size_in_bytes_hpp__
