//
// Copyright (c) 2026 INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/serialization.hpp"
#endif // PINOCCHIO_LSP

namespace boost
{
  namespace serialization
  {

    namespace internal
    {
      template<typename Scalar, int Options, std::size_t Alignment>
      struct BlockDiagonalMatrixAccessor
      : public ::pinocchio::BlockDiagonalMatrixTpl<Scalar, Options, Alignment>
      {
        typedef ::pinocchio::BlockDiagonalMatrixTpl<Scalar, Options, Alignment> Base;
        using Base::m_cols;
        using Base::m_matrix_block_elements;
        using Base::m_matrix_stack;
        using Base::m_rows;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options, std::size_t Alignment>
    void serialize(
      Archive & ar,
      ::pinocchio::BlockDiagonalMatrixTpl<Scalar, Options, Alignment> & matrix,
      const unsigned int /*version*/)
    {
      typedef internal::BlockDiagonalMatrixAccessor<Scalar, Options, Alignment> Accessor;
      Accessor & matrix_ = static_cast<Accessor &>(matrix);
      ar & make_nvp("m_rows", matrix_.m_rows);
      ar & make_nvp("m_cols", matrix_.m_cols);
      ar & make_nvp("m_matrix_stack", matrix_.m_matrix_stack);

      auto & m_matrix_block_elements = matrix_.m_matrix_block_elements;
      ar & make_nvp("m_matrix_block_elements", m_matrix_block_elements);
      std::vector<std::size_t> index_map;
      index_map.reserve(m_matrix_block_elements.size());
      if (Archive::is_loading::value)
      {
        ar & make_nvp("index_map", index_map);
        for (size_t k = 0; k < index_map.size(); ++k)
        {
          const auto index_value = index_map[k];
          if (index_value == std::numeric_limits<std::size_t>::max())
            continue;
          auto & block = m_matrix_block_elements[k];

          block.remap(matrix_.m_matrix_stack[index_value]);
        }
      }
      else
      {
        const auto & m_matrix_stack = matrix_.m_matrix_stack;

        for (const auto & block : m_matrix_block_elements)
        {
          const auto block_data = block.data();
          if (block_data != nullptr)
          {
            const auto it = std::find_if(
              m_matrix_stack.begin(), m_matrix_stack.end(),
              [&block_data](const auto & stack_elt) { return block_data == stack_elt.data(); });

            assert(it != m_matrix_stack.end() && "must_never happened");

            std::size_t stack_elt_index = std::size_t(std::distance(m_matrix_stack.begin(), it));
            index_map.push_back(stack_elt_index);
          }
          else
            index_map.push_back(std::numeric_limits<std::size_t>::max());
        }

        ar & make_nvp("index_map", index_map);
      }
    }

  } // namespace serialization
} // namespace boost