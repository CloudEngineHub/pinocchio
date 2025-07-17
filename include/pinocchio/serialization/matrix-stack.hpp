//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_matrix_stack_hpp__
#define __pinocchio_serialization_matrix_stack_hpp__

#include "pinocchio/serialization/eigen.hpp"

#include "pinocchio/container/matrix-stack.hpp"

namespace boost
{
  namespace serialization
  {

    namespace internal
    {
      template<typename MatrixLike, std::size_t Alignment>
      struct MatrixStackAccessor : public ::pinocchio::MatrixStackTpl<MatrixLike, Alignment>
      {
        typedef ::pinocchio::MatrixStackTpl<MatrixLike, Alignment> Base;
        using Base::m_data_ptr;
        using Base::m_matrix_maps;
        using Base::m_memory_capacity;
        using Base::m_offsets;

        using Base::free;
        using Base::incr_ptr;
        using Base::malloc;
      };
    } // namespace internal

    template<typename Archive, typename MatrixLike, std::size_t Alignment>
    void serialize(
      Archive & ar,
      ::pinocchio::MatrixStackTpl<MatrixLike, Alignment> & matrix_stack,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::MatrixStackTpl<MatrixLike, Alignment> MatrixStack;
      typedef typename MatrixStack::MapType MapType;
      typedef typename MatrixStack::Scalar Scalar;
      typedef internal::MatrixStackAccessor<MatrixLike, Alignment> Accessor;
      Accessor & matrix_stack_ = static_cast<Accessor &>(matrix_stack);

      auto & offsets = matrix_stack_.m_offsets;
      auto & memory_capacity = matrix_stack_.m_memory_capacity;
      // auto & matrix_maps = matrix_stack_.m_matrix_maps;

      if (Archive::is_loading::value)
        matrix_stack.clear();

      std::size_t stack_size = matrix_stack.size();
      ar & make_nvp("stack_size", stack_size);
      if (Archive::is_loading::value)
      {
        for (std::size_t i = 0; i < stack_size; ++i)
        {
          // const auto offset = offsets[i];
          Eigen::Index rows = -1, cols = -1;
          ar & make_nvp("rows", rows);
          ar & make_nvp("cols", cols);
          matrix_stack.push_back(rows, cols);
          auto & matrix_map = matrix_stack.back();
          ar & make_nvp("map", matrix_map);
        }
      }
      else // writting mode
      {
        for (std::size_t i = 0; i < matrix_stack.size(); ++i)
        {
          auto & matrix_map = matrix_stack[i];
          Eigen::Index rows = matrix_map.rows(), cols = matrix_map.cols();
          ar & make_nvp("rows", rows);
          ar & make_nvp("cols", cols);

          ar & make_nvp("map", matrix_map);
        }
      }
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_matrix_stack_hpp__
