//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_zero_cone_hpp__
#define __pinocchio_algorithm_constraints_zero_cone_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename _Scalar, int _Options>
  struct traits<ZeroConeTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options
    };
    typedef FullSpaceConeTpl<Scalar, _Options> DualCone;
  };

  ///  \brief Null set containing (0 singleton).
  template<typename _Scalar, int _Options>
  struct ZeroConeTpl : ConeBase<ZeroConeTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef ConeBase<ZeroConeTpl> Base;
    typedef typename traits<ZeroConeTpl>::DualCone DualCone;

    /// \brief Cast to base class.
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    /// \brief Const cast to base class.
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    /// \brief Comparison operator
    bool operator==(const ZeroConeTpl & other) const
    {
      return base() == other.base();
    }

    /// \brief Difference  operator
    bool operator!=(const ZeroConeTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Returns the dual cone of this.
    DualCone dual() const
    {
      return DualCone();
    }

    using Base::isInside;
    /// \brief Check whether a vector x is zero.
    ///
    /// \param[in] f vector to check (assimilated to a  force vector).
    ///
    template<typename VectorLike>
    bool isInsideImpl(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      return x.isZero(prec);
    }

    using Base::project;
    /// \brief Project a vector x into set.
    ///
    /// \param[in] x a vector to project.
    /// \param[in] res result of the projection.
    ///
    template<typename VectorLikeIn, typename VectorLikeOut>
    void projectImpl(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeOut> & res_) const
    {
      PINOCCHIO_UNUSED_VARIABLE(x);
      auto & res = res_.const_cast_derived();
      res.setZero();
    }

  }; // ZeroConeTpl

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_zero_cone_hpp__
