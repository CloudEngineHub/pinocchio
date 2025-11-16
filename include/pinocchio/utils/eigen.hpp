//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_eigen_hpp__
#define __pinocchio_utils_eigen_hpp__

#include "pinocchio/utils/fwd.hpp"
#include <algorithm>

namespace pinocchio
{

  namespace helper
  {
    template<class T>
    struct is_eigen_noalias : std::false_type
    {
    };

    template<typename ExpressionType, template<typename> class StorageBase>
    struct is_eigen_noalias<Eigen::NoAlias<ExpressionType, StorageBase>> : std::true_type
    {
    };

    template<class T>
    inline constexpr bool is_eigen_noalias_v =
      is_eigen_noalias<std::remove_cv_t<std::remove_reference_t<T>>>::value;

  } // namespace helper
  namespace internal
  {

    template<typename ExpressionType, template<typename> class StorageBase>
    class PromoteStaticOp
    {
    public:
      typedef typename ExpressionType::Scalar Scalar;

      EIGEN_DEVICE_FUNC PromoteStaticOp(ExpressionType & expression)
      : m_expression(expression)
      {
      }
      EIGEN_DEVICE_FUNC PromoteStaticOp(ExpressionType && expression)
      : m_expression(expression)
      {
      }

      template<typename OtherDerived>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType &
      operator=(const StorageBase<OtherDerived> & other)
      {
        call_assignment_no_alias(
          m_expression, other.derived(),
          Eigen::internal::assign_op<Scalar, typename OtherDerived::Scalar>());
        return m_expression;
      }

      template<typename Lhs, typename Rhs, int Option>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType &
      operator=(const Eigen::Product<Lhs, Rhs, Option> & matrix_product)
      {
        const auto & lhs = matrix_product.lhs();
        const auto & rhs = matrix_product.rhs();

        typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Lhs) PlainLhs;
        typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Rhs) PlainRhs;
        const auto lhs_map = make_map<PlainLhs>(lhs);
        const auto rhs_map = make_map<PlainRhs>(rhs);
        const auto matrix_map_product = lhs_map * rhs_map;
        typedef decltype(matrix_map_product) OtherDerived;

        if constexpr (helper::is_eigen_noalias_v<ExpressionType>)
        {
          typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(
            std::remove_cv_t<std::remove_reference_t<decltype(expression().expression())>>)
            PlainExpression;
          auto result_matrix_map = make_map<PlainExpression>(expression().expression());

          // result_matrix_map.noalias() = matrix_map_product;
          Eigen::internal::call_assignment_no_alias(
            result_matrix_map, matrix_map_product,
            Eigen::internal::assign_op<
              typename PlainExpression::Scalar, typename OtherDerived::Scalar>());
        }
        else
        {
          typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(ExpressionType) PlainExpression;
          auto result_matrix_map = make_map<PlainExpression>(expression());

          // result_matrix_map = matrix_map_product;
          Eigen::internal::call_assignment(
            result_matrix_map, matrix_map_product,
            Eigen::internal::assign_op<
              typename PlainExpression::Scalar, typename OtherDerived::Scalar>());
        }

        return m_expression;
      }
      // template <typename OtherDerived>
      // EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType& operator+=(const
      // StorageBase<OtherDerived>& other) {
      //   call_assignment_no_alias(m_expression, other.derived(),
      //                            internal::add_assign_op<Scalar, typename
      //                            OtherDerived::Scalar>());
      //   return m_expression;
      // }

      // template <typename OtherDerived>
      // EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType& operator-=(const
      // StorageBase<OtherDerived>& other) {
      //   call_assignment_no_alias(m_expression, other.derived(),
      //                            internal::sub_assign_op<Scalar, typename
      //                            OtherDerived::Scalar>());
      //   return m_expression;
      // }

      EIGEN_DEVICE_FUNC ExpressionType & expression() const
      {
        return m_expression;
      }

    protected:
      ExpressionType & m_expression;

      typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> DynamicStride;

      template<typename MatrixLike, typename MatrixDerived>
      static Eigen::Map<MatrixLike, 0, DynamicStride>
      make_map(const Eigen::MatrixBase<MatrixDerived> & _plain_object)
      {
        auto & plain_object = _plain_object.const_cast_derived();
        const DynamicStride stride = {plain_object.outerStride(), plain_object.innerStride()};
        return {plain_object.data(), plain_object.rows(), plain_object.cols(), stride};
      }
    };

  } // namespace internal

  template<typename MatrixExpression>
  internal::PromoteStaticOp<MatrixExpression, Eigen::MatrixBase>
  promote_static_op(const Eigen::MatrixBase<MatrixExpression> & matrix_expression)
  {
    return {matrix_expression.const_cast_derived()};
  }

  template<typename MatrixExpression, template<typename> class StorageBase>
  internal::PromoteStaticOp<Eigen::NoAlias<MatrixExpression, StorageBase>, Eigen::MatrixBase>
  promote_static_op(Eigen::NoAlias<MatrixExpression, StorageBase> && matrix_expression)
  {
    return {std::forward<Eigen::NoAlias<MatrixExpression, StorageBase>>(matrix_expression)};
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_utils_eigen_hpp__
