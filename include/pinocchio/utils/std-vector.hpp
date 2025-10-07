//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_std_vector_hpp__
#define __pinocchio_utils_std_vector_hpp__

#include <vector>

#include "pinocchio/utils/template-template-parameter.hpp"

namespace pinocchio
{
  namespace internal
  {
    template<typename V>
    struct std_vector_extract_allocator_type
    {
      template<typename T>
      using type =
        typename extract_template_template_parameter<typename V::allocator_type>::template type<T>;
    };

    template<typename V>
    struct std_vector_with_same_allocator
    {
      template<typename T>
      using allocator_type = typename std_vector_extract_allocator_type<V>::template type<T>;

      template<typename T>
      using type = std::vector<T, allocator_type<T>>;
    };

  } // namespace internal

  template<typename T, typename Allocator, class Func>
  void apply_for_each(std::vector<T, Allocator> & vector, const Func & func)
  {
    std::for_each(vector.begin(), vector.end(), func);
  }

  /**
   * @brief Creates a vector of holder objects that wrap the elements of a given vector.
   *
   * This function takes a reference to a `std::vector` of elements of type `T`
   * and constructs a new `std::vector` containing `Holder<T>` objects,
   * each created from the corresponding element in the input vector.
   *
   * Typical use case: producing a vector of `std::reference_wrapper<T>` or
   * other holder objects for easy element access or reference semantics.
   *
   * @tparam T         The element type stored in the input vector.
   * @tparam Allocator The allocator type used by the input vector.
   * @tparam Holder    A class template that accepts a single template parameter
   *                   (for example, `std::reference_wrapper` or a custom holder template).
   * @param vec        Reference to the vector containing elements of type `T`.
   *
   * @return A new vector of type `std::vector<Holder<T>>`, where each element
   *         wraps or references the corresponding element from the input vector.
   *
   * @note If `Holder` is a reference wrapper (e.g. `std::reference_wrapper`),
   *       the returned holders will refer to the original elements in `vec`.
   *       Make sure the lifetime of `vec` exceeds the lifetime of the returned vector
   *       to avoid dangling references.
   *
   * @see std::reference_wrapper
   */
  template<typename T, typename Allocator, template<typename H> class Holder>
  std::vector<Holder<T>> make_held_vector(std::vector<T, Allocator> & vec)
  {
    typedef std::vector<Holder<T>> WrappedTVector;
    return WrappedTVector(vec.cbegin(), vec.cend());
  }

  /**
   * @brief Creates a vector of holder objects that wrap the elements of a given const vector.
   *
   * This function takes a constant reference to a `std::vector` of elements of type `T`
   * and returns a new `std::vector` containing `Holder<const T>` objects constructed
   * from each element of the input vector.
   *
   * Typical use case: producing a vector of `std::reference_wrapper<const T>` or
   * another lightweight holder type from a vector of const elements.
   *
   * @tparam T         The element type stored (const-qualified) in the input vector.
   * @tparam Allocator The allocator type used by the input vector.
   * @tparam Holder    A class template that accepts a single type parameter
   *                   (e.g., `std::reference_wrapper` or a custom holder template).
   * @param vec        The input vector containing elements of type `const T`.
   *
   * @return A new vector of type `std::vector<Holder<const T>>`, where each element
   *         wraps the corresponding element from the input vector.
   *
   * @note The elements are copied or wrapped using the constructor of `Holder<const T>`
   *       that takes a `const T&`. To avoid dangling references, ensure the lifetime
   *       of the original elements outlives the returned holders if `Holder` is a reference
   * wrapper.
   *
   * @see std::reference_wrapper
   */
  template<typename T, typename Allocator, template<typename H> class Holder>
  std::vector<Holder<const T>> make_held_vector(const std::vector<const T, Allocator> & vec)
  {
    typedef std::vector<Holder<const T>> WrappedTVector;
    return WrappedTVector(vec.cbegin(), vec.cend());
  }
} // namespace pinocchio

#endif // __pinocchio_utils_std_vector_hpp__
