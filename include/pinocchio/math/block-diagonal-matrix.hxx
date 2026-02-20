//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_math_block_diagonal_matrix_hxx__
#define __pinocchio_math_block_diagonal_matrix_hxx__

#include "pinocchio/math/matrix-info.hpp"
#include "pinocchio/math/matrix.hpp"
#include "pinocchio/alloca.hpp"

namespace pinocchio
{

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVectorType>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::init_or_rebuild(
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression)
  {
    // typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(DiagonalVectorType) PlainDiagonalVectorType;
    // const auto diagonal_terms = make_map<ConstMatrixMap>(diagonal_expression.diagonal());
    // ConstMatrixBlockElement block_info = {
    //   pinocchio::MatrixBlockType::Diagonal, diagonal_terms.size(), diagonal_terms};
    const auto & diagonal_terms = diagonal_expression.diagonal();
    ConstMatrixBlockElement block_info = {
      pinocchio::MatrixBlockType::Diagonal, diagonal_terms.size()};

    const std::vector<ConstMatrixBlockElement> input_block_pattern = {{block_info}};

    init_or_rebuild(input_block_pattern);
    m_matrix_block_elements.back().container() = diagonal_terms;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVectorType>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::BlockDiagonalMatrixTpl(
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression)
  {
    init_or_rebuild(diagonal_expression);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVectorType>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::rebuild(
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression)
  {
    init_or_rebuild(diagonal_expression);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::BlockDiagonalMatrixTpl(
    const std::vector<MatrixBlockElement> & input_block_pattern)
  {
    init_or_rebuild(input_block_pattern);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename _MatrixBlockElement>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::rebuild(
    const std::vector<_MatrixBlockElement> & new_block_pattern)
  {
    init_or_rebuild(new_block_pattern);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename _MatrixBlockElement>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::rebuild(
    const _MatrixBlockElement * new_block_pattern, const size_t size)
  {
    init_or_rebuild(new_block_pattern, size);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::clear()
  {
    m_rows = m_cols = -1;
    m_matrix_stack.clear();
    m_matrix_block_elements.clear();
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename _MatrixBlockElement>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::init_or_rebuild(
    const std::vector<_MatrixBlockElement> & input_block_pattern)
  {
    init_or_rebuild(input_block_pattern.data(), input_block_pattern.size());
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename _MatrixBlockElement>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::init_or_rebuild(
    const _MatrixBlockElement * input_block_pattern, const size_t size)
  {
    clear();

    static_assert(
      pinocchio::internal::is_specialization_of_v<
        _MatrixBlockElement, pinocchio::MatrixBlockElementTpl>,
      "_MatrixBlockElement is not of type pinocchio::MatrixBlockElementTpl<...>");

    // analysis block pattern and extract memory/size info
    const std::size_t num_blocks = size;
    MatrixInfo * memory_block_sizes =
      static_cast<MatrixInfo *>(PINOCCHIO_ALLOCA(num_blocks * sizeof(MatrixInfo)));
    std::size_t memory_block_id = 0;

    m_matrix_block_elements.reserve(size);
    m_rows = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
      const auto & block_info = input_block_pattern[i];
      assert(block_info.type() != MatrixBlockType::Undefined);

      m_rows += block_info.size();

      m_matrix_block_elements.push_back({block_info.type(), block_info.size()});

      if (block_info.type() == MatrixBlockType::ScalarIdentity)
      {
        const MatrixInfo memory_block_size = {1, 1};
        memory_block_sizes[memory_block_id++] = memory_block_size;
      }
      else if (block_info.type() == MatrixBlockType::Diagonal)
      {
        const MatrixInfo memory_block_size = {block_info.size(), 1};
        memory_block_sizes[memory_block_id++] = memory_block_size;
      }
      else if (block_info.type() == MatrixBlockType::Plain)
      {
        const MatrixInfo memory_block_size = {block_info.size(), block_info.size()};
        memory_block_sizes[memory_block_id++] = memory_block_size;
      }
    }

    m_matrix_stack.rebuild(memory_block_sizes, memory_block_id);

    // Fill with data
    std::size_t matrix_stack_id = 0;
    for (std::size_t i = 0; i < m_matrix_block_elements.size(); ++i)
    {
      auto & block_info = m_matrix_block_elements[i];
      auto & input_block_info = input_block_pattern[i];

      assert(block_info.type() != MatrixBlockType::Undefined);
      if (isDataBlock(block_info.type()))
      {
        auto matrix_map = m_matrix_stack[matrix_stack_id];
        block_info.remap(matrix_map); // remap data to the matrix stack
        if (input_block_info.data() != nullptr)
          block_info.container() = input_block_info.container(); // copy data
        // otherwise, unitialized
        matrix_stack_id++;
      }
    }

    m_cols = m_rows;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename MatrixDerived>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::evalTo(
    const Eigen::MatrixBase<MatrixDerived> & _matrix) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _matrix.rows(), rows(), "The input matrix has not the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _matrix.cols(), cols(), "The input matrix has not the right number of columns.");

    _matrix.const_cast_derived().setZero();
    assign_op<internal::assign_op>(_matrix.const_cast_derived());
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename MatrixDerived>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::addTo(
    const Eigen::MatrixBase<MatrixDerived> & _matrix) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _matrix.rows(), rows(), "The input matrix has not the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _matrix.cols(), cols(), "The input matrix has not the right number of columns.");

    assign_op<internal::add_assign_op>(_matrix.const_cast_derived());
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename MatrixDerived>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::subTo(
    const Eigen::MatrixBase<MatrixDerived> & _matrix) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _matrix.rows(), rows(), "The input matrix has not the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _matrix.cols(), cols(), "The input matrix has not the right number of columns.");

    assign_op<internal::sub_assign_op>(_matrix.const_cast_derived());
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename AssignOp, typename MatrixDerived>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::assign_op(
    const Eigen::MatrixBase<MatrixDerived> & _matrix) const
  {
    auto & matrix = _matrix.const_cast_derived();

    Eigen::Index row_id = 0;
    for (const auto & matrix_block_elt : m_matrix_block_elements)
    {
      const auto size = matrix_block_elt.size();
      auto matrix_block = matrix.block(row_id, row_id, size, size);

      matrix_block_elt.template assign_op<AssignOp>(matrix_block);
      row_id += size;
    }

    assert(row_id == cols());
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename MatrixDerived>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::matrix(
    const Eigen::MatrixBase<MatrixDerived> & _matrix) const
  {
    evalTo(_matrix.const_cast_derived());
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  typename BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::Matrix
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::matrix() const
  {
    Matrix res(rows(), cols());
    matrix(res);
    return res;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  bool BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::hasNaN() const
  {
    for (const auto & block_info : m_matrix_block_elements)
    {
      if (block_info.hasNaN())
        return true;
    }
    return false;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename AssignOp, typename MatrixDerivedRhs, typename MatrixDerivedRes>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::applyOnTheRight(
    const Eigen::MatrixBase<MatrixDerivedRhs> & rhs,
    const Eigen::MatrixBase<MatrixDerivedRes> & _res) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      rhs.rows(), cols(), "The input rhs matrix has not the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      rhs.cols(), _res.cols(), "The input rhs and res matrices has not the right number of cols.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _res.rows(), rows(), "The input res matrix has not the right number of rows.");

    auto & res = _res.const_cast_derived();

    Eigen::Index row_id = 0;
    for (const auto & block_info : m_matrix_block_elements)
    {
      const auto block_size = block_info.size();

      const auto rhs_block = rhs.middleRows(row_id, block_size);
      auto res_block = res.middleRows(row_id, block_size);

      assert(
        rhs_block.data() != res_block.data()
        && "There is aliasing between rhs_block and res_block.");

      switch (block_info.type())
      {
      case MatrixBlockType::Zero: {
        AssignOp::run(Matrix::Zero(rhs_block.rows(), rhs_block.cols()), res_block);
        break;
      }
      case MatrixBlockType::Identity: {
        AssignOp::run(rhs_block, res_block);
        break;
      }
      case MatrixBlockType::ScalarIdentity: {
        const auto & map = block_info.map;
        const auto & scalar = map(0, 0);
        AssignOp::run(scalar * rhs_block, res_block);
        break;
      }
      case MatrixBlockType::Diagonal: {
        const auto & map = block_info.map;
        AssignOp::run(map.asDiagonal() * rhs_block, res_block.noalias());
        break;
      }
      case MatrixBlockType::Plain: {
        const auto & map = block_info.map;
        AssignOp::run(map * rhs_block, res_block.noalias());
        break;
      }
      default:
        assert(false && "Should never happen");
      }

      row_id += block_size;
    }
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename AssignOp, typename MatrixDerivedLhs, typename MatrixDerivedRes>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::applyOnTheLeft(
    const Eigen::MatrixBase<MatrixDerivedLhs> & lhs,
    const Eigen::MatrixBase<MatrixDerivedRes> & _res) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      lhs.cols(), rows(), "The input lhs matrix has not the right number of cols.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      lhs.rows(), _res.rows(), "The input lhs and res matrices has not the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      _res.cols(), cols(), "The input res matrix has not the right number of cols.");

    auto & res = _res.const_cast_derived();

    Eigen::Index col_id = 0;
    for (const auto & block_info : m_matrix_block_elements)
    {
      const auto block_size = block_info.size();

      const auto lhs_block = lhs.middleCols(col_id, block_size);
      auto res_block = res.middleCols(col_id, block_size);

      assert(
        lhs_block.data() != res_block.data()
        && "There is aliasing between lhs_block and res_block.");

      switch (block_info.type())
      {
      case MatrixBlockType::Zero: {
        AssignOp::run(Matrix::Zero(lhs_block.rows(), lhs_block.cols()), res_block);
        break;
      }
      case MatrixBlockType::Identity: {
        AssignOp::run(lhs_block, res_block);
        break;
      }
      case MatrixBlockType::ScalarIdentity: {
        const auto & map = block_info.map;
        const auto & scalar = map(0, 0);
        AssignOp::run(lhs_block * scalar, res_block);
        break;
      }
      case MatrixBlockType::Diagonal: {
        const auto & map = block_info.map;
        AssignOp::run(lhs_block * map.asDiagonal(), res_block.noalias());
        break;
      }
      case MatrixBlockType::Plain: {
        const auto & map = block_info.map;
        AssignOp::run(lhs_block * map, res_block.noalias());
        break;
      }
      default:
        assert(false && "Should never happen");
      }

      col_id += block_size;
    }
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  Inverse<BlockDiagonalMatrixTpl<Scalar, Options, Alignment>>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::inverse() const
  {
    return {*this};
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment> &
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::operator=(
    const BlockDiagonalMatrixTpl & other)
  {
    if (this == &other)
      return *this;

    m_rows = other.m_rows;
    m_cols = other.m_cols;
    m_matrix_stack = other.m_matrix_stack;
    m_matrix_block_elements = other.m_matrix_block_elements;

    size_t matrix_stack_id = 0;
    for (auto & block_info : m_matrix_block_elements)
    {
      if (isDataBlock(block_info.type()))
      {
        auto matrix_map = m_matrix_stack[matrix_stack_id];
        block_info.remap(matrix_map);
        matrix_stack_id++;
      }
    }

    return *this;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVectorType>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment> &
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::operator=(
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression)
  {
    rebuild(diagonal_expression);
    return *this;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVectorType>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment> &
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::operator+=(
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression)
  {
    Sum<
      BlockDiagonalMatrixTpl<Scalar, Options, Alignment>,
      Eigen::DiagonalWrapper<DiagonalVectorType>>(*this, diagonal_expression)
      .evalTo(*this);
    return *this;
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVectorType>
  Sum<
    BlockDiagonalMatrixTpl<Scalar, Options, Alignment>,
    Eigen::DiagonalWrapper<DiagonalVectorType>>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::operator+(
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression) const
  {
    return {*this, diagonal_expression};
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::Zero(const Eigen::Index size)
  {
    MatrixBlockElement block_info = {pinocchio::MatrixBlockType::Zero, size};
    const std::vector<MatrixBlockElement> input_block_pattern = {{block_info}};

    return BlockDiagonalMatrixTpl(input_block_pattern);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>
  BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::ScalarIdentity(
    const Eigen::Index size, const Scalar & value)
  {
    typedef Eigen::Matrix<Scalar, 1, 1> M11;
    M11 value_mat = M11(value);
    const auto matrix_map = make_map<MatrixMap>(value_mat);
    MatrixBlockElement block_info = {pinocchio::MatrixBlockType::ScalarIdentity, size, matrix_map};
    const std::vector<MatrixBlockElement> input_block_pattern = {{block_info}};

    return BlockDiagonalMatrixTpl(input_block_pattern);
  }

  template<typename Scalar, int Options, std::size_t Alignment>
  template<typename DiagonalVector>
  void BlockDiagonalMatrixTpl<Scalar, Options, Alignment>::diagonal(
    const Eigen::MatrixBase<DiagonalVector> & _diagonal_elements) const
  {
    auto & diagonal_elements = _diagonal_elements.const_cast_derived();
    Eigen::Index row_id = 0;
    for (const auto & block_info : m_matrix_block_elements)
    {
      const auto block_size = block_info.size();
      auto diagonal_elements_segment = diagonal_elements.segment(row_id, block_size);
      block_info.diagonal(diagonal_elements_segment);
      row_id += block_size;
    }
  }

} // namespace pinocchio

#endif // #ifndef __pinocchio_math_block_diagonal_matrix_hxx__
