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
    /**
     * @brief Helper trait for reusing the allocator type of a given `std::vector`‑like type.
     *
     * This trait extracts the allocator template template‑parameter from a type `V`
     * (typically a `std::vector<...>`) and exposes a nested alias template `type`
     * that can be used to rebind the allocator to a different element type.
     *
     * Internally it relies on another trait,
     * `extract_template_template_parameter`, to retrieve the allocator’s
     * template template‑parameter so that it can be applied as
     * `Allocator<T>` for an arbitrary element type `T`.
     *
     * Example:
     * @code
     *   using VecDouble = std::vector<double, std::allocator<double>>;
     *   using Extractor = std_vector_extract_allocator_type<VecDouble>;
     *
     *   // Get the allocator type suitable for 'int'
     *   using AllocForInt = Extractor::type<int>;
     *   std::vector<int, AllocForInt> intVec;
     * @endcode
     *
     * @tparam V The vector‑like type from which to extract the allocator.
     *
     * @see extract_template_template_parameter
     */
    template<typename V>
    struct std_vector_extract_allocator_type
    {
      template<typename T>
      using type =
        typename extract_template_template_parameter<typename V::allocator_type>::template type<T>;
    };

    /**
     * @brief Helper trait for creating `std::vector` types that reuse the allocator of another
     * vector.
     *
     * This metafunction extracts the allocator type from an existing vector-like type `V`
     * and defines nested aliases to build new `std::vector` specializations using the same
     * allocator but possibly with a different element type.
     *
     * The typical purpose is to ensure allocator consistency when generating vectors of
     * related element types (e.g., converting a `std::vector<double>` into a
     * `std::vector<int>` that uses the same allocator configuration).
     *
     * Example:
     * @code
     *   using VecDouble = std::vector<double, std::allocator<double>>;
     *   using VectorAllocatorAdapter = std_vector_with_same_allocator<VecDouble>;
     *
     *   // Defines a std::vector<int, std::allocator<int>>
     *   using VecInt = VectorAllocatorAdapter::type<int>;
     *
     *   // Alternatively, the allocator alone:
     *   using AllocForInt = VectorAllocatorAdapter::allocator_type<int>;
     * @endcode
     *
     * @tparam V  The existing vector-like type whose allocator type should be reused.
     *
     * @see std_vector_extract_allocator_type
     */
    template<typename V>
    struct std_vector_with_same_allocator
    {
      template<typename T>
      using allocator_type = typename std_vector_extract_allocator_type<V>::template type<T>;

      template<typename T>
      using type = std::vector<T, allocator_type<T>>;
    };

  } // namespace internal

  /**
   * @brief Applies a given function to each element in a std::vector.
   *
   * This function uses `std::for_each` to apply the provided function
   * to each element in the input vector.
   *
   * @tparam T The type of elements stored in the vector.
   * @tparam Allocator The allocator used by the vector.
   * @tparam Func The type of the function to be applied.
   *
   * @param vector The vector whose elements the function will be applied to.
   * @param func The function to apply to each element. It should accept a single argument of type
   * `T&`.
   */
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
   * @tparam Holder    A class template that accepts a single template parameter
   *                   (for example, `std::reference_wrapper` or a custom holder template).
   * @tparam T         The element type stored in the input vector.
   * @tparam Allocator The allocator type used by the input vector.
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
  template<template<typename H> class Holder, typename T, typename Allocator>
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
   * @tparam Holder    A class template that accepts a single type parameter
   *                   (e.g., `std::reference_wrapper` or a custom holder template).
   * @tparam T         The element type stored (const-qualified) in the input vector.
   * @tparam Allocator The allocator type used by the input vector.
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
  template<template<typename H> class Holder, typename T, typename Allocator>
  std::vector<Holder<const T>> make_held_vector(const std::vector<T, Allocator> & vec)
  {
    typedef std::vector<Holder<const T>> WrappedTVector;
    return WrappedTVector(vec.cbegin(), vec.cend());
  }
} // namespace pinocchio

#endif // __pinocchio_utils_std_vector_hpp__
