//
// Copyright (c) 2026 INRIA
//
#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include <type_traits>
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  template<typename T>
  struct EigenMatrixExpression;

  template<typename Op, typename Type>
  struct UnaryOperator;

  template<typename Op, typename LhsType, typename RhsType>
  struct BinaryOperator;

  template<typename T>
  struct Inverse;

  template<typename Lhs, typename Rhs>
  struct Sum;

  template<typename Scalar, int Options, std::size_t Alignment = alignof(std::max_align_t)>
  struct BlockDiagonalMatrixTpl;

  template<typename MatrixOrMap, typename Enable = void>
  struct MatrixBlockElementTpl;

  template<typename Derived>
  struct MatrixBlockElementPlain;

  template<typename T>
  struct is_floating_point : ::std::is_floating_point<T>
  {
  };

  ///  \brief Foward declaration of TaylorSeriesExpansion.
  template<typename Scalar>
  struct TaylorSeriesExpansion;

} // namespace pinocchio
