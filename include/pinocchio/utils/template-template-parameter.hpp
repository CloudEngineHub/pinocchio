//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_template_template_parameter_hpp__
#define __pinocchio_utils_template_template_parameter_hpp__

namespace pinocchio
{
  namespace internal
  {
    /**
     * @brief Primary template for extracting a template template parameter from a class template.
     *
     * This struct serves as a general fallback for types that are not recognized as instantiations
     * of class templates. It simply aliases the original type `C` as a nested template.
     *
     * @tparam C A class type (possibly a class template instantiation or not).
     *
     * This template is typically specialized for cases where `C` is a template instantiation.
     */
    template<typename C>
    struct extract_template_template_parameter
    {
      /**
       * @brief A nested alias that defaults to the original type `C`.
       *
       * This template alias acts as a passthrough. It exists for completeness and fallback
       * behavior.
       *
       * @tparam Unused Dummy parameter for SFINAE compatibility.
       */
      template<typename = void>
      using type = C;
    };

    /**
     * @brief Specialization of extract_template_template_parameter for class template
     * instantiations.
     *
     * This specialization extracts the template template parameter `C` from a class template
     * instantiation of the form `C<Parameters...>`. It exposes a nested template alias `type`
     * that can be used to re-instantiate `C` with different parameters.
     *
     * @tparam C The template template parameter (i.e., a class template taking variadic
     * parameters).
     * @tparam Parameters The parameter pack used in the original instantiation of `C`.
     *
     * @code
     * template<template<typename...> class Tpl, typename... Args>
     * using Extracted = typename extract_template_template_parameter<Tpl<Args...>>::template
     * type<NewArgs...>;
     * @endcode
     */
    template<template<class...> class C, class... Parameters>
    struct extract_template_template_parameter<C<Parameters...>>
    {
      /**
       * @brief A nested alias to re-instantiate the extracted template `C` with different
       * parameters.
       *
       * This allows applying the extracted template template parameter `C` to a new set of template
       * arguments.
       *
       * @tparam Other The new set of template parameters to instantiate `C` with.
       */
      template<class... Other>
      using type = C<Other...>;
    };

  } // namespace internal
} // namespace pinocchio

#endif // __pinocchio_utils_template_template_parameter_hpp__
