//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_sets_box_set_hpp__
#define __pinocchio_algorithm_constraints_sets_box_set_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/sets/set-base.hpp"
#include "pinocchio/math/matrix.hpp"

namespace pinocchio
{

  template<typename _Scalar, int _Options>
  struct traits<BoxSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
  };

  /// \brief Box set defined by a lower and an upper bounds [lb;ub].
  ///  This operator does not own any data.
  ///  Instead, it points to the data owned by some other struct.
  ///  Creating and copying an instance of this struct is free.
  template<typename _Scalar, int _Options>
  struct BoxSetTpl : SetBase<BoxSetTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef SetBase<BoxSetTpl> Base;

    using Base::project;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    /// \brief Constructor lb and ub.
    /// Internally, a const reference to these vectors is kept.
    ///
    /// \param[in] lb box lower bound
    /// \param[in] ub box upper bound
    BoxSetTpl(const Vector & lb, const Vector & ub)
    : lb(lb)
    , ub(ub)
    {
      assert(lb.size() == ub.size());
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        arrayCompareAll(lb, ub, internal::ComparisonOperators::LE),
        "Some components of lb are greater than ub");
    }

    /// \brief Generic constructor which can take any `Parameter` struct.
    /// Creates a link between the `mu` in `Parameter` and this->mu.
    ///
    /// \param[in] params Generic parameters, must contain the field `mu`.
    template<typename Parameters>
    BoxSetTpl(const Parameters & params)
    : BoxSetTpl(params.lb, params.ub)
    {
    }

    /// \brief Copy constructor.
    BoxSetTpl(const BoxSetTpl & other) = default;

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

    /// \brief Copy operator
    BoxSetTpl & operator=(const BoxSetTpl & other) = default;

    /// \brief Comparison operator
    bool operator==(const BoxSetTpl & other) const
    {
      return base() == other.base() && lb == other.lb && ub == other.ub;
    }

    /// \brief Difference operator
    bool operator!=(const BoxSetTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Check whether lb <= ub for all components
    bool isValid() const
    {
      return (lb.array() <= ub.array()).all();
    }

    /// \brief Project the value given as input for the given row index.
    Scalar rowiseProject(const Eigen::Index row_id, const Scalar value) const
    {
      assert(row_id < lb.size());
      return math::max(lb[row_id], math::min(ub[row_id], value));
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc Base::isInside
    template<typename VectorLike>
    bool isInsideImpl(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
    {
      assert(check_expression_if_real<Scalar>(prec >= 0) && "prec should be positive");
      return check_expression_if_real<Scalar>((x - project(x)).norm() <= prec);
    }

    /// \copydoc Base::project
    template<typename VectorLikeIn, typename VectorLikeOut>
    void projectImpl(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeOut> & res_) const
    {
      pinocchio::arrayBound(x, lb, ub, res_);
    }

    /// \copydoc Base::scaledProject
    template<typename VectorLikeIn, typename VectorLikeIn2, typename VectorLikeOut>
    void scaledProjectImpl(
      const Eigen::MatrixBase<VectorLikeIn> & x,
      const Eigen::MatrixBase<VectorLikeIn2> & scale,
      const Eigen::MatrixBase<VectorLikeOut> & res_) const
    {
      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
      assert((scale.array() > 0).all() && "scale vector should be positive");
      pinocchio::arrayBound(
        x, (lb.array() / scale.array()).matrix(), (ub.array() / scale.array()).matrix(), res_);
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    }

    // ------------------------------
    // MEMBERS
    // ------------------------------

    /// \brief Reference to the lower bound of the box.
    const Vector & lb;
    /// \brief Reference to the upper bound of the box.
    const Vector & ub;

  }; // BoxSetTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_sets_box_set_hpp__
