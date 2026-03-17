//
// Copyright (c) 2026 INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/math.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<typename Scalar, int Options, std::size_t Alignment, typename DiagonalVectorType>
  struct Sum<
    BlockDiagonalMatrixTpl<Scalar, Options, Alignment>,
    Eigen::DiagonalWrapper<DiagonalVectorType>>
  : BlockDiagonalMatrixExpression<Sum<
      BlockDiagonalMatrixTpl<Scalar, Options, Alignment>,
      Eigen::DiagonalWrapper<DiagonalVectorType>>>
  {
    typedef BlockDiagonalMatrixTpl<Scalar, Options, Alignment> LhsType;
    typedef Eigen::DiagonalWrapper<DiagonalVectorType> RhsType;

    Sum(const LhsType & lhs, const RhsType & rhs)
    : m_lhs(lhs)
    , m_rhs(rhs)
    {
    }

    /// @brief \copydoc Base::evalTo
    template<typename ResScalar, int ResOptions, std::size_t ResAlignment>
    void evalTo(BlockDiagonalMatrixTpl<ResScalar, ResOptions, ResAlignment> & res) const
    {
      typedef BlockDiagonalMatrixTpl<ResScalar, ResOptions, ResAlignment> ResType;

      const auto & diag = rhs().diagonal();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(
        diag.size(), lhs().rows(),
        "The size of the diagonal expression does not match the number of rows of the block "
        "diagonal matrix.");

      // Check if rebuild is needed
      bool need_rebuild = false;
      if (
        res.rows() != lhs().rows() || res.cols() != lhs().cols()
        || res.blocks().size() != lhs().blocks().size())
      {
        need_rebuild = true;
      }
      else
      {
        for (std::size_t i = 0; i < lhs().blocks().size(); ++i)
        {
          // Upgrading rule: anything not Diagonal/Plain becomes Diagonal
          const auto current_type = lhs().blocks()[i].type();
          const auto target_type =
            (current_type != MatrixBlockType::Diagonal && current_type != MatrixBlockType::Plain)
              ? MatrixBlockType::Diagonal
              : current_type;

          if (res.blocks()[i].type() != target_type)
          {
            need_rebuild = true;
            break;
          }
        }
      }

      ResType * res_ptr = &res;
      ResType tmp_res;

      if (need_rebuild)
      {
        typedef Eigen::Matrix<ResScalar, Eigen::Dynamic, Eigen::Dynamic, ResOptions> ResMatrix;
        typedef Eigen::Map<ResMatrix, to_eigen_alignment(ResAlignment)> ResMatrixMap;
        typedef MatrixBlockElementTpl<ResMatrixMap> ResMatrixBlockElement;

        typedef typename ResType::MatrixBlockElement MatrixBlockElement;
        static_assert(
          pinocchio::internal::is_same_type<MatrixBlockElement, ResMatrixBlockElement>::value,
          "MatrixBlockElement is not of type pinocchio::MatrixBlockElementTpl<Eigen::Map>");

        const std::size_t num_blocks = lhs().blocks().size();
        MatrixBlockElement * new_pattern = static_cast<MatrixBlockElement *>(
          PINOCCHIO_ALLOCA(num_blocks * sizeof(MatrixBlockElement)));

        for (std::size_t i = 0; i < num_blocks; ++i)
        {
          const auto & block = lhs().blocks()[i];
          MatrixBlockType new_type = block.type();
          if (new_type != MatrixBlockType::Diagonal && new_type != MatrixBlockType::Plain)
            new_type = MatrixBlockType::Diagonal;
          new (new_pattern + i) MatrixBlockElement(new_type, block.size());
        }

        if (&res == &lhs())
        {
          tmp_res.rebuild(new_pattern, num_blocks);
          res_ptr = &tmp_res;
        }
        else
        {
          res.rebuild(new_pattern, num_blocks);
        }
      }

      Eigen::Index row_id = 0;
      for (std::size_t i = 0; i < lhs().blocks().size(); ++i)
      {
        const auto & lhs_block = lhs().blocks()[i];
        auto & res_block = res_ptr->blocks()[i];
        const auto block_size = lhs_block.size();
        const auto diag_segment = diag.segment(row_id, block_size);

        if (!need_rebuild && &res == &lhs())
        {
          // In-place: res += rhs
          if (res_block.type() == MatrixBlockType::Diagonal)
          {
            res_block.container() += diag_segment;
          }
          else
          {
            res_block.container().diagonal() += diag_segment;
          }
        }
        else
        {
          // res = lhs + rhs
          if (res_block.type() == MatrixBlockType::Diagonal)
          {
            assert(lhs_block.type() != MatrixBlockType::Plain);
            lhs_block.diagonal(res_block.container());
            res_block.container() += diag_segment;
          }
          else
          {
            assert(lhs_block.type() == MatrixBlockType::Plain);
            lhs_block.evalTo(res_block.container());
            res_block.container().diagonal() += diag_segment;
          }
        }

        row_id += block_size;
      }

      if (need_rebuild && &res == &lhs())
      {
        res = std::move(tmp_res);
      }
    }

    const LhsType & lhs() const
    {
      return m_lhs;
    }

    const RhsType & rhs() const
    {
      return m_rhs;
    }

    Eigen::Index rows() const
    {
      return m_lhs.rows();
    }

    Eigen::Index cols() const
    {
      return m_lhs.cols();
    }

    Eigen::Index size() const
    {
      return m_lhs.size();
    }

  protected:
    const LhsType & m_lhs; // block diagonal matrix ref
    const RhsType m_rhs;   // diagonal expression
  };

} // namespace pinocchio
