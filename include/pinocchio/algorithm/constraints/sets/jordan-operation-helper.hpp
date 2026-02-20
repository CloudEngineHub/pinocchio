//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_algorithm_constraints_sets_jordan_operations_helper_hpp__
#define __pinocchio_algorithm_constraints_sets_jordan_operations_helper_hpp__

#include <type_traits>

namespace pinocchio
{
  namespace helper
  {
    // --- does traits<T>::JordanOperation exist?
    template<class T, class = void>
    struct has_jordan_operation : std::false_type
    {
    };

    template<class T>
    struct has_jordan_operation<T, std::void_t<typename traits<T>::JordanOperation>>
    : std::true_type
    {
    };

    // --- check that the jordan operation is on the non negative orthant
    template<class T, class = void>
    struct has_non_negative_orthant_jordan_operation : std::false_type
    {
    };

    template<class T>
    struct has_non_negative_orthant_jordan_operation<
      T,
      std::void_t<typename traits<T>::JordanOperation>>
    : std::is_same<
        typename traits<T>::JordanOperation,
        NonNegativeOrthantJordanOperationTpl<typename traits<T>::Scalar, traits<T>::Options>>
    {
    };

    // --- check that the jordan operation is on the second order cone
    template<class T, class = void>
    struct has_second_order_cone_jordan_operation : std::false_type
    {
    };

    template<class T>
    struct has_second_order_cone_jordan_operation<
      T,
      std::void_t<typename traits<T>::JordanOperation>>
    : std::is_same<
        typename traits<T>::JordanOperation,
        SecondOrderConeJordanOperationTpl<typename traits<T>::Scalar, traits<T>::Options>>
    {
    };

    // --- check that the jordan operation is on the zero symmetric cone
    template<class T, class = void>
    struct has_zero_cone_jordan_operation : std::false_type
    {
    };

    template<class T>
    struct has_zero_cone_jordan_operation<T, std::void_t<typename traits<T>::JordanOperation>>
    : std::is_same<
        typename traits<T>::JordanOperation,
        ZeroConeJordanOperationTpl<typename traits<T>::Scalar, traits<T>::Options>>
    {
    };
  } // namespace helper

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_sets_jordan_operations_helper_hpp__
