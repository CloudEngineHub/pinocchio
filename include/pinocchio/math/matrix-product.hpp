//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_math_matrix_product_hpp__
#define __pinocchio_math_matrix_product_hpp__

#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{

  template<template<typename, typename> class EigenOp, typename Lhs, typename Rhs, typename Res>
  void matrix_product(
    const Eigen::MatrixBase<Lhs> & lhs,
    const Eigen::MatrixBase<Rhs> & rhs,
    const Eigen::MatrixBase<Res> & res)
  {
    using Scalar = typename Rhs::Scalar;

    const auto rows = lhs.rows();
    const auto cols = rhs.cols();
    const auto inner_dim = lhs.cols();

    const auto * lhs_data = lhs.derived().data();
    const auto * rhs_data = rhs.derived().data();
    auto * res_data = res.const_cast_derived().data();

    auto lhs_index = [&](const Eigen::Index i, const Eigen::Index j) -> Eigen::Index {
      if constexpr (Lhs::IsRowMajor)
        return j + i * inner_dim;
      else // lhs(i,k) = lhs_data[k*rows + i]
        return j * rows + i;
    };

    auto rhs_index = [&](const Eigen::Index i, const Eigen::Index j) -> Eigen::Index {
      if constexpr (Rhs::IsRowMajor)
        return j + i * cols;
      else // rhs(k,j) = rhs_data[j*inner_dim + k]
        return j * inner_dim + i;
    };

    auto res_index = [&](const Eigen::Index i, const Eigen::Index j) -> Eigen::Index {
      if constexpr (Res::IsRowMajor)
        return j + i * cols;
      else // res(i,j) = res_data[j*rows + i]
        return j * rows + i;
    };

#pragma omp simd
    for (Eigen::Index j = 0; j < cols; ++j) // loop over columns of result
    {
      for (Eigen::Index i = 0; i < rows; ++i) // loop over rows of result
      {
        Scalar sum = Scalar(0);
        for (Eigen::Index k = 0; k < inner_dim; ++k)
        {
          sum += lhs_data[lhs_index(i, k)] * rhs_data[rhs_index(k, j)];
        }
        typedef EigenOp<typename Res::Scalar, typename Lhs::Scalar> Op;
        if constexpr (internal::is_specialization_of_v<Op, Eigen::internal::assign_op>)
        {
          res_data[res_index(i, j)] = sum;
        }
        else if constexpr (internal::is_specialization_of_v<Op, Eigen::internal::add_assign_op>)
        {
          res_data[res_index(i, j)] += sum;
        }
        else if constexpr (internal::is_specialization_of_v<Op, Eigen::internal::sub_assign_op>)
        {
          res_data[res_index(i, j)] -= sum;
        }
      }
    }
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_math_matrix_product_hpp__
