//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_size_in_bytes_hpp__
#define __pinocchio_utils_size_in_bytes_hpp__

#include <cstddef>
#include <type_traits>
#include <utility> // for std::declval

#include "pinocchio/utils/fwd.hpp"

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

  /**
   * @brief Helper struct providing a fallback implementation to compute the size (in bytes)
   *        of a given object or type.
   *
   * Specialize this struct for custom types that do not have a @c sizeInBytes() member function,
   * in order to define how their size should be computed.
   *
   * @tparam T Type of the object whose size in bytes will be computed.
   */
  template<typename T>
  struct sizeInBytesImpl
  {
    /**
     * @brief Compute the size in bytes of a given object value of type @p T.
     *
     * This static method should return the number of bytes occupied by @p value.
     * The default implementation may rely on sizeof(T), or it may be specialized
     * for custom data structures that store variable-length data.
     *
     * @param value The object whose size in bytes should be computed.
     * @return The size in bytes of @p value.
     */
    static std::size_t run(const T & value);
  };

  /**
   * @brief Compute the size in bytes of an object @p value.
   *
   * This function first checks whether the type @p T provides a member function
   * named @c sizeInBytes(). If it does, that method is called.
   * Otherwise, it falls back to calling the static implementation
   * provided by the @c sizeInBytesImpl<T> struct.
   *
   * Typical use case:
   * @code
   * MyType obj;
   * std::size_t sz = sizeInBytes(obj);
   * @endcode
   *
   * @tparam T Type of the input object.
   * @param value The input object whose size in bytes will be computed.
   * @return The size in bytes of @p value.
   */
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

  /**
   * @brief Get the size in bytes of a type @p T at compile time.
   *
   * This overload simply returns @c sizeof(T), which is sufficient for
   * trivial or fixed-size types.
   *
   * Typical use case:
   * @code
   * std::size_t sz = sizeInBytes<int>(); // returns sizeof(int)
   * @endcode
   *
   * @tparam T Type whose compile-time size in bytes is requested.
   * @return The size in bytes of type @p T.
   */
  template<typename T>
  std::size_t sizeInBytes()
  {
    return sizeof(T);
  }

} // namespace pinocchio

#endif // __pinocchio_utils_size_in_bytes_hpp__
