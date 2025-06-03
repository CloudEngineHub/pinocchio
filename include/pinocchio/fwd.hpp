//
// Copyright (c) 2018-2024 CNRS INRIA
//

#ifndef __pinocchio_fwd_hpp__
#define __pinocchio_fwd_hpp__

// Forward declaration of the main pinocchio namespace
namespace pinocchio
{
}

#ifdef _WIN32
  #include <windows.h>
  #undef far
  #undef near
#endif

#include <cassert>

#ifdef PINOCCHIO_EIGEN_CHECK_MALLOC
  #ifndef EIGEN_RUNTIME_NO_MALLOC
    #define EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
    #define EIGEN_RUNTIME_NO_MALLOC
  #endif
#endif

#include "pinocchio/macros.hpp"
#include "pinocchio/deprecated.hpp"
#include "pinocchio/warning.hpp"
#include "pinocchio/config.hpp"
#include "pinocchio/unsupported.hpp"

// Include Eigen components
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
  #include <Eigen/AccelerateSupport>
#endif

#include "pinocchio/eigen-macros.hpp"
#ifdef PINOCCHIO_WITH_EIGEN_TENSOR_MODULE
  #include <unsupported/Eigen/CXX11/Tensor>
#endif

#include "pinocchio/utils/helpers.hpp"
#include "pinocchio/utils/cast.hpp"
#include "pinocchio/utils/check.hpp"

#include "pinocchio/container/boost-container-limits.hpp"

#ifdef PINOCCHIO_EIGEN_CHECK_MALLOC
  #ifndef EIGEN_RUNTIME_NO_MALLOC
    #define EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
    #define EIGEN_RUNTIME_NO_MALLOC
  #endif
#endif

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
  #include <Eigen/AccelerateSupport>
#endif

#ifdef PINOCCHIO_EIGEN_CHECK_MALLOC
  #ifdef EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
    #undef EIGEN_RUNTIME_NO_MALLOC
    #undef EIGEN_RUNTIME_NO_MALLOC_WAS_NOT_DEFINED
  #endif
#endif

#include "pinocchio/core/binary-op.hpp"
#include "pinocchio/core/unary-op.hpp"

#include <cstddef> // std::size_t

namespace pinocchio
{
  ///
  /// \brief Common traits structure to fully define base classes for CRTP.
  ///
  template<class C>
  struct traits
  {
  };

  /// \brief Blank type
  struct Blank
  {
  };

  namespace internal
  {
    template<typename T>
    struct traits
    {
    };
  } // namespace internal

  template<class Derived>
  struct NumericalBase
  {
    typedef typename traits<Derived>::Scalar Scalar;
  };

  ///
  /// \brief Type of the cast of a class C templated by Scalar and Options, to a new NewScalar type.
  ///        This class should be specialized for each types.
  ///
  template<typename NewScalar, class C>
  struct CastType;

  ///
  ///  \brief Cast scalar type from type FROM to type TO.
  ///
  template<typename To, typename From>
  struct ScalarCast
  {
    static To cast(const From & value)
    {
      return static_cast<To>(value);
    }
  };

  template<typename To, typename From>
  To scalar_cast(const From & value)
  {
    return ScalarCast<To, From>::cast(value);
  }

  /// \brief Argument position.
  ///        Used as template parameter to refer to an argument.
  enum ArgumentPosition
  {
    ARG0 = 0,
    ARG1 = 1,
    ARG2 = 2,
    ARG3 = 3,
    ARG4 = 4
  };

  /// \brief Assignment operator list.
  ///
  enum AssignmentOperatorType
  {
    SETTO,
    ADDTO,
    RMTO
  };

  ///  \brief Assignment operator tags
  template<AssignmentOperatorType val>
  struct AssignmentOperatorTag
  {
  };

  using SetTo = AssignmentOperatorTag<SETTO>;
  using AddTo = AssignmentOperatorTag<ADDTO>;
  using RmTo = AssignmentOperatorTag<RMTO>;

  /** This value means that a positive quantity (e.g., a size) is not known at compile-time, and
   * that instead the value is stored in some runtime variable.
   */
  const int Dynamic = -1;

  /// \brief Undefined return type
  ///        This is an helper structure to help internal diagnosis.
  struct UndefinedReturnType;

  typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;

  namespace internal
  {
    /**
     * @brief Type trait to determine whether a type is a specialization of a given class template.
     *
     * This trait evaluates to `std::true_type` if the first argument is an instantiation of
     * the template provided as the second argument. Otherwise, it evaluates to `std::false_type`.
     *
     * @tparam T The type to test (e.g., `std::vector<int>`).
     * @tparam Template The class template to test against (e.g., `std::vector`).
     *
     * ### Example
     * @code
     * static_assert(is_specialization_of<std::vector<int>, std::vector>::value, "Should be true");
     * static_assert(!is_specialization_of<int, std::vector>::value, "Should be false");
     * @endcode
     *
     * This works for any template that can be expressed as `template<typename...> class`.
     */
    template<typename T, template<typename...> class Template>
    struct is_specialization_of : std::false_type
    {
    };

    /**
     * @brief Partial specialization: case when @p T is an instantiation of @p Template.
     *
     * This specialization activates when the first template parameter is of the form
     * `Template<Args...>`, making the trait derive from `std::true_type`.
     *
     * @tparam Template The class template to test against.
     * @tparam Args Template argument pack used in the specialization being tested.
     */
    template<template<typename...> class Template, typename... Args>
    struct is_specialization_of<Template<Args...>, Template> : std::true_type
    {
    };

    /**
     * @brief Helper variable template for @ref is_specialization_of.
     *
     * Provides a `constexpr bool` that simplifies access to the result of
     * `is_specialization_of<T, Template>::value`.
     *
     * This is analogous to standard library variable templates like
     * `std::is_same_v` or `std::is_base_of_v`.
     *
     * @tparam T The type to test (e.g., `std::vector<int>`).
     * @tparam Template The class template to test against (e.g., `std::vector`).
     *
     * ### Example
     * @code
     * static_assert(is_specialization_of_v<std::vector<int>, std::vector>);
     * static_assert(!is_specialization_of_v<int, std::vector>);
     * @endcode
     *
     * @see is_specialization_of
     */
    template<typename T, template<typename...> class Template>
    inline constexpr bool is_specialization_of_v = is_specialization_of<T, Template>::value;
  } // namespace internal
} // namespace pinocchio

#include "pinocchio/context.hpp"
#include "pinocchio/alloca.hpp"

#endif // #ifndef __pinocchio_fwd_hpp__
