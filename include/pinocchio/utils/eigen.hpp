//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_utils_eigen_hpp__
#define __pinocchio_utils_eigen_hpp__

#include "pinocchio/utils/fwd.hpp"
#include <type_traits>

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

    template<class T>
    struct is_eigen_product : std::false_type
    {
    };

    template<typename Lhs, typename Rhs, int Option>
    struct is_eigen_product<Eigen::Product<Lhs, Rhs, Option>> : std::true_type
    {
    };

    template<class T>
    inline constexpr bool is_eigen_product_v =
      is_eigen_product<std::remove_cv_t<std::remove_reference_t<T>>>::value;

    template<typename T>
    struct remove_eigen_noalias
    {
      typedef T type;
    };

    template<typename ExpressionType, template<typename> class StorageBase>
    struct remove_eigen_noalias<Eigen::NoAlias<ExpressionType, StorageBase>>
    {
      typedef ExpressionType type;
    };

    template<typename T, typename = void>
    inline constexpr bool has_fixed_rows_v = false;

    template<typename T>
    inline constexpr bool has_fixed_rows_v<T, std::void_t<decltype(T::RowsAtCompileTime)>> =
      (T::RowsAtCompileTime != Eigen::Dynamic);

    template<typename T, typename = void>
    inline constexpr bool has_fixed_cols_v = false;

    template<typename T>
    inline constexpr bool has_fixed_cols_v<T, std::void_t<decltype(T::ColsAtCompileTime)>> =
      (T::ColsAtCompileTime != Eigen::Dynamic);

    template<typename T, typename = void>
    inline constexpr bool has_fixed_size_v = false;

    template<typename T>
    inline constexpr bool has_fixed_size_v<
      T,
      std::void_t<decltype(T::RowsAtCompileTime), decltype(T::ColsAtCompileTime)>> =
      has_fixed_rows_v<T> && has_fixed_cols_v<T>;

  } // namespace helper

  namespace internal
  {

    template<typename ExpressionType, template<typename> class StorageBase>
    class PromoteStaticOp
    {
    public:
      typedef typename ExpressionType::Scalar Scalar;

      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(
        typename helper::remove_eigen_noalias<ExpressionType>::type) PlainExpression;

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
        return dispatch<Eigen::internal::assign_op>(matrix_product.derived());
      }

      template<typename Lhs, typename Rhs, int Option>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType &
      operator+=(const Eigen::Product<Lhs, Rhs, Option> & matrix_product)
      {
        return dispatch<Eigen::internal::add_assign_op>(matrix_product.derived());
      }

      template<typename Lhs, typename Rhs, int Option>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType &
      operator-=(const Eigen::Product<Lhs, Rhs, Option> & matrix_product)
      {
        return dispatch<Eigen::internal::sub_assign_op>(matrix_product.derived());
      }

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

      template<typename MatrixLike, typename MatrixDerived, template<typename> class _StorageBase>
      static Eigen::Map<MatrixLike, 0, DynamicStride>
      make_map(const Eigen::NoAlias<MatrixDerived, _StorageBase> & _plain_object_noalias)
      {
        return make_map<MatrixLike>(_plain_object_noalias.expression().const_cast_derived());
      }

      template<
        template<typename, typename> class Op,
        typename Dst,
        template<typename> class _StorageBase,
        typename Src>
      inline constexpr void
      call_assignment(Eigen::NoAlias<Dst, _StorageBase> && dst, const Src & src)
      {
        typedef typename Dst::Scalar S1;
        typedef typename Src::Scalar S2;

        Eigen::internal::call_assignment_no_alias(dst.expression(), src, Op<S1, S2>());
      }

      template<template<typename, typename> class Op, typename Dst, typename Src>
      inline constexpr void call_assignment(const Eigen::MatrixBase<Dst> & dst, const Src & src)
      {
        typedef typename Dst::Scalar S1;
        typedef typename Src::Scalar S2;

        Eigen::internal::call_assignment(dst.const_cast_derived(), src, Op<S1, S2>());
      }

      template<typename Result, typename Lhs, typename Rhs>
      struct MatrixProductDimensions
      {
        static constexpr int RowsAtCompileTime = Lhs::RowsAtCompileTime != Eigen::Dynamic
                                                   ? static_cast<int>(Lhs::RowsAtCompileTime)
                                                   : static_cast<int>(Result::RowsAtCompileTime);
        static constexpr int ColsAtCompileTime = Rhs::ColsAtCompileTime != Eigen::Dynamic
                                                   ? static_cast<int>(Rhs::ColsAtCompileTime)
                                                   : static_cast<int>(Result::ColsAtCompileTime);
        static constexpr int InnerDimensionAtCompileTime =
          Lhs::ColsAtCompileTime != Eigen::Dynamic ? static_cast<int>(Lhs::ColsAtCompileTime)
                                                   : static_cast<int>(Rhs::RowsAtCompileTime);
      };

      template<typename Product>
      static constexpr bool is_static_size_product()
      {
        static_assert(helper::is_eigen_product_v<Product>);
        using Lhs = typename Product::Lhs;
        using Rhs = typename Product::Rhs;

        typedef MatrixProductDimensions<PlainExpression, Lhs, Rhs> Dims;
        constexpr int RowsAtCompileTime = Dims::RowsAtCompileTime;
        constexpr int ColsAtCompileTime = Dims::ColsAtCompileTime;
        constexpr int InnerDimensionAtCompileTime = Dims::InnerDimensionAtCompileTime;

        return RowsAtCompileTime != Eigen::Dynamic && ColsAtCompileTime != Eigen::Dynamic
               && InnerDimensionAtCompileTime != Eigen::Dynamic;
      }

      template<typename Product>
      static constexpr bool is_partial_static_size_product()
      {
        static_assert(helper::is_eigen_product_v<Product>);
        using Lhs = typename Product::Lhs;
        using Rhs = typename Product::Rhs;

        typedef MatrixProductDimensions<PlainExpression, Lhs, Rhs> Dims;
        constexpr int RowsAtCompileTime = Dims::RowsAtCompileTime;
        constexpr int ColsAtCompileTime = Dims::ColsAtCompileTime;
        constexpr int InnerDimensionAtCompileTime = Dims::InnerDimensionAtCompileTime;

        return (RowsAtCompileTime != Eigen::Dynamic && ColsAtCompileTime != Eigen::Dynamic)
               || (ColsAtCompileTime != Eigen::Dynamic && InnerDimensionAtCompileTime != Eigen::Dynamic)
               || (InnerDimensionAtCompileTime != Eigen::Dynamic && RowsAtCompileTime != Eigen::Dynamic);
      }

      template<template<typename, typename> class EigenOp, typename Lhs, typename Rhs, int Option>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ExpressionType &
      dispatch(const Eigen::Product<Lhs, Rhs, Option> & matrix_product)
      {
        typedef Eigen::Product<Lhs, Rhs, Option> ProductType;
        if constexpr (is_static_size_product<ProductType>())
        {
          static_dispatch<EigenOp>(matrix_product);
        }
        else
        {
          dynamic_dispatch<EigenOp>(matrix_product);
        }
        return m_expression;
      }

      template<template<typename, typename> class EigenOp, typename Lhs, typename Rhs, int Option>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
      dynamic_dispatch(const Eigen::Product<Lhs, Rhs, Option> & matrix_product)
      {
        static_dispatch<EigenOp>(matrix_product);
      }

      template<template<typename, typename> class EigenOp, typename Lhs, typename Rhs, int Option>
      EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
      static_dispatch(const Eigen::Product<Lhs, Rhs, Option> & matrix_product)
      {
        const auto & lhs = matrix_product.lhs();
        const auto & rhs = matrix_product.rhs();

        typedef MatrixProductDimensions<PlainExpression, Lhs, Rhs> Dims;
        constexpr int RowsAtCompileTime = Dims::RowsAtCompileTime;
        constexpr int ColsAtCompileTime = Dims::ColsAtCompileTime;
        constexpr int InnerDimensionAtCompileTime = Dims::InnerDimensionAtCompileTime;

        typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Lhs) PlainLhs;
        typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Rhs) PlainRhs;

        typedef Eigen::Matrix<
          typename PlainLhs::Scalar, RowsAtCompileTime, InnerDimensionAtCompileTime,
          PlainLhs::Options>
          PlainLhsStaticSize;
        typedef Eigen::Matrix<
          typename PlainRhs::Scalar, InnerDimensionAtCompileTime, ColsAtCompileTime,
          PlainRhs::Options>
          PlainRhsStaticSize;

        const auto lhs_map = make_map<PlainLhsStaticSize>(lhs);
        const auto rhs_map = make_map<PlainRhsStaticSize>(rhs);
        const auto matrix_map_product = lhs_map * rhs_map;
        // typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE((Eigen::Product<Lhs, Rhs, Option>))
        // PlainProductResult;

        // std::cout << "PlainProductResult: " << typeid(PlainProductResult).name() << std::endl;
        // std::cout << "RowsAtCompileTime: " << RowsAtCompileTime << std::endl;
        // std::cout << "ColsAtCompileTime: " << ColsAtCompileTime << std::endl;
        // std::cout << "InnerDimensionAtCompileTime: " << InnerDimensionAtCompileTime << std::endl;

        typedef Eigen::Matrix<
          typename PlainExpression::Scalar, RowsAtCompileTime, ColsAtCompileTime,
          PlainExpression::Options>
          PlainExpressionStaticSize;
        auto result_matrix_map = make_map<PlainExpressionStaticSize>(expression());

        if constexpr (helper::is_eigen_noalias_v<ExpressionType>)
        {
          call_assignment<EigenOp>(result_matrix_map.noalias(), matrix_map_product);
        }
        else
        {
          call_assignment<EigenOp>(result_matrix_map, matrix_map_product);
        }
      }

    }; // struct PromoteStaticOp

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
