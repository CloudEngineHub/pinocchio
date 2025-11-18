//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_orthant_cone_hpp__
#define __pinocchio_algorithm_constraints_orthant_cone_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename _Scalar>
  struct traits<NonNegativeOrthantConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef NonNegativeOrthantConeTpl<Scalar> DualCone;
  };

  template<typename _Scalar>
  struct NonNegativeOrthantConeTpl : ConeBase<NonNegativeOrthantConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef ConeBase<NonNegativeOrthantConeTpl> Base;
    typedef typename traits<NonNegativeOrthantConeTpl>::DualCone DualCone;

    /// \brief Cast to base class.
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    /// \brief Const cast to base class.
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    using Base::derived;
    /// \brief Returns the dual cone associated with this.
    ///
    /// \remarks Orthant cone are by definition self dual.
    DualCone dual() const
    {
      return derived();
    }

    /// \brief Comparison operator
    bool operator==(const NonNegativeOrthantConeTpl & other) const
    {
      return base() == other.base();
    }

    /// \brief Difference  operator
    bool operator!=(const NonNegativeOrthantConeTpl & other) const
    {
      return !(*this == other);
    }

    using Base::isInside;
    /// \brief Check whether a vector x lies within the orthant.
    ///
    /// \param[in] x vector to check .
    ///
    template<typename VectorLike>
    bool isInsideImpl(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(prec >= 0 && "prec should be positive");
      return (x - project(x)).norm() <= prec;
    }

    using Base::project;
    /// \brief Project a vector x into orthant.
    ///
    /// \param[in] x a vector to project.
    /// \param[in] res result of the projection.
    ///
    template<typename VectorLikeIn, typename VectorLikeOut>
    void projectImpl(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeOut> & res_) const
    {
      res_.const_cast_derived() = x.array().max(Scalar(0)).matrix();
    }

    /// \brief Project the value given as input for the given row index.
    Scalar rowiseProject(const Eigen::DenseIndex /*row_id*/, const Scalar value) const
    {
      return math::max(Scalar(0), value);
    }

  }; // struct NonNegativeOrthantConeTpl

} // namespace pinocchio
#endif // ifndef __pinocchio_algorithm_constraints_orthant_cone_hpp__
