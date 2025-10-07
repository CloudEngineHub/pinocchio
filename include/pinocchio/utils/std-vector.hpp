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

  template<typename T, typename Allocator, template<typename H> class Holder>
  std::vector<Holder<T>> make_held_vector(std::vector<T, Allocator> & vec)
  {
    typedef std::vector<Holder<T>> WrappedTVector;
    return WrappedTVector(vec.cbegin(), vec.cend());
  }

  template<typename T, typename Allocator, template<typename H> class Holder>
  std::vector<Holder<const T>> make_held_vector(const std::vector<const T, Allocator> & vec)
  {
    typedef std::vector<Holder<const T>> WrappedTVector;
    return WrappedTVector(vec.cbegin(), vec.cend());
  }
} // namespace pinocchio

#endif // __pinocchio_utils_std_vector_hpp__
