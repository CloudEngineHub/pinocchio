//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_serialization_matrix_block_element_hpp__
#define __pinocchio_serialization_matrix_block_element_hpp__

#include "pinocchio/serialization/eigen.hpp"
#include "pinocchio/math/matrix-block-element.hpp"

namespace boost
{
  namespace serialization
  {

    namespace internal
    {
      template<typename Matrix>
      struct MatrixBlockElementTplAccessor : public ::pinocchio::MatrixBlockElementTpl<Matrix>
      {
        typedef ::pinocchio::MatrixBlockElementTpl<Matrix> Base;
        using Base::m_size;
        using Base::m_type;
      };
    } // namespace internal

    template<typename Archive, typename Matrix>
    void serialize(
      Archive & ar,
      ::pinocchio::MatrixBlockElementTpl<Matrix> & _matrix_block_element,
      const unsigned int /*version*/)
    {
      typedef internal::MatrixBlockElementTplAccessor<Matrix> Accessor;
      Accessor & matrix_block_element = static_cast<Accessor &>(_matrix_block_element);

      ar & make_nvp("type", matrix_block_element.m_type);
      ar & make_nvp("size", matrix_block_element.m_size);

      if constexpr (pinocchio::helper::is_eigen_matrix_v<Matrix>)
      {
        auto & container = matrix_block_element.container();
        ar & make_nvp("container", container);
      }
    }

  } // namespace serialization
} // namespace boost

#endif // ifndef __pinocchio_serialization_matrix_block_element_hpp__
