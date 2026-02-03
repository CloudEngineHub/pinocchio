//
// Copyright (c) 2025 INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/math.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  namespace internal
  {
    template<>
    struct MatrixInversionCodeGeneratedImpl<3>
    {
      template<typename M1, typename M2>
      static EIGEN_STRONG_INLINE void
      run(const Eigen::MatrixBase<M1> & matrix, const Eigen::MatrixBase<M2> & matrix_inverse_)
      {
        typedef typename M1::Scalar Scalar;

        assert(is_symmetric(matrix, math::sqrt(dummy_precision<Scalar>())));

        const auto & input_vec = matrix.reshaped();
        auto & matrix_inverse = matrix_inverse_.const_cast_derived();
        auto output_vec = matrix_inverse.reshaped();
      }
    };
  } // namespace internal
} // namespace pinocchio