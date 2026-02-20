//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_math_block_diagonal_matrix_expression_hpp__
#define __pinocchio_math_block_diagonal_matrix_expression_hpp__

namespace pinocchio
{
  template<typename Derived>
  struct BlockDiagonalMatrixExpression : BlockDiagonalMatrixBase<Derived>
  {
    typedef BlockDiagonalMatrixBase<Derived> Base;
    using Base::derived;

    /// @brief Evaluates this expression and stores it in res.
    template<typename Scalar, int Options, std::size_t Alignment>
    void evalTo(BlockDiagonalMatrixTpl<Scalar, Options, Alignment> & res) const
    {
      derived().evalTo(res.derived());
    }
  }; // struct BlockDiagonalMatrixExpression

} // namespace pinocchio

#endif // #ifndef __pinocchio_math_block_diagonal_matrix_expression_hpp__
