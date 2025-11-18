//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_full_space_cone_hpp__
#define __pinocchio_algorithm_constraints_full_space_cone_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename _Scalar, int _Options>
  struct traits<FullSpaceConeTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef ZeroConeTpl<Scalar, _Options> DualCone;
  };

  ///  \brief Unbounded set covering the whole space
  template<typename _Scalar, int _Options>
  struct FullSpaceConeTpl : ConeBase<FullSpaceConeTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef ConeBase<FullSpaceConeTpl> Base;
    typedef typename traits<FullSpaceConeTpl>::DualCone DualCone;

    /// \brief Cast operator
    template<typename NewScalar>
    FullSpaceConeTpl<NewScalar, Options> cast() const
    {
      typedef FullSpaceConeTpl<NewScalar, Options> ReturnType;
      return ReturnType();
    }

    /// \brief Cast to base class
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    /// \brief Const cast to base class
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    /// \brief Returns the dual cone of this.
    DualCone dual() const
    {
      return DualCone();
    }

    /// \brief Comparison operator
    bool operator==(const FullSpaceConeTpl & other) const
    {
      return base() == other.base();
    }

    /// \brief Difference  operator
    bool operator!=(const FullSpaceConeTpl & other) const
    {
      return !(*this == other);
    }

    using Base::isInside;
    /// \brief Check whether a vector x lies within the set.
    /// Any vector x always belong the the unbounded set.
    ///
    /// \param[in] f vector to check (assimilated to a  force vector).
    ///
    template<typename VectorLike>
    bool isInsideImpl(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      PINOCCHIO_UNUSED_VARIABLE(x);
      PINOCCHIO_UNUSED_VARIABLE(prec);
      return true;
    }

    using Base::project;
    /// \brief Project a vector x into the set.
    ///
    /// \param[in] x a vector to project.
    /// \param[in] res result of the projection.
    ///
    template<typename VectorLikeIn, typename VectorLikeOut>
    void projectImpl(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeOut> & res_) const
    {
      res_.const_cast_derived() = x;
    }

  }; // FullSpaceConeTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_full_space_cone_hpp__
