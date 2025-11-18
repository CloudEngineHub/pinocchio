//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_zero_cone_hpp__
#define __pinocchio_algorithm_constraints_zero_cone_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, ZeroConeTpl<Scalar, Options>>
  {
    typedef ZeroConeTpl<NewScalar, Options> type;
  };

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

    /// \brief Constructor from a given size
    ///
    explicit ZeroConeTpl(const Eigen::DenseIndex size)
    : m_size(size)
    {
    }

    /// \brief Copy constructor.
    ZeroConeTpl(const ZeroConeTpl & other) = default;

    /// \brief Copy operator
    ZeroConeTpl & operator=(const ZeroConeTpl & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    ZeroConeTpl<NewScalar, Options> cast() const
    {
      typedef ZeroConeTpl<NewScalar, Options> ReturnType;
      return ReturnType(this->size());
    }

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
      return base() == other.base() && m_size == other.m_size;
    }

    /// \brief Difference  operator
    bool operator!=(const ZeroConeTpl & other) const
    {
      return !(*this == other);
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

    /// \brief Returns the dimension of the ambiant space.
    Eigen::DenseIndex dim() const
    {
      return m_size;
    }

    Eigen::DenseIndex size() const
    {
      return m_size;
    }

    /// \brief Resize by calling the resize method of Eigen.
    void resize(Eigen::DenseIndex new_size)
    {
      m_size = new_size;
    }

    /// \brief Resize by calling the conservativeResize method of Eigen.
    void conservativeResize(Eigen::DenseIndex new_size)
    {
      this->resize(new_size);
    }

    DualCone dual() const
    {
      return DualCone(m_size);
    }

  protected:
    Eigen::DenseIndex m_size;
  }; // ZeroConeTpl

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_zero_cone_hpp__
