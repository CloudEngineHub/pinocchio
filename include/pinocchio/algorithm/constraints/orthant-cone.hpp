//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_orthant_cone_hpp__
#define __pinocchio_algorithm_constraints_orthant_cone_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/algorithm/constraints/cone-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar>
  struct CastType<NewScalar, OrthantConeTpl<Scalar>>
  {
    typedef OrthantConeTpl<NewScalar> type;
  };

  template<typename _Scalar>
  struct traits<OrthantConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef OrthantConeTpl<Scalar> DualCone;
  };

  template<typename _Scalar>
  struct OrthantConeTpl : ConeBase<OrthantConeTpl<_Scalar>>
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef ConeBase<OrthantConeTpl> Base;
    typedef typename traits<OrthantConeTpl>::DualCone DualCone;

    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    /// \brief Default constructor
    ///
    OrthantConeTpl()
    : m_size(0)
    {
    }

    /// \brief Constructor from a given size
    ///
    explicit OrthantConeTpl(const Eigen::DenseIndex size)
    : m_size(size)
    {
    }

    /// \brief Copy constructor.
    OrthantConeTpl(const OrthantConeTpl & other) = default;

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, OrthantConeTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, OrthantConeTpl>::type ReturnType;
      return ReturnType(size());
    }

    /// \brief Copy operator
    OrthantConeTpl & operator=(const OrthantConeTpl & other) = default;

    /// \brief Comparison operator
    bool operator==(const OrthantConeTpl & other) const
    {
      return base() == other.base() && m_size == other.m_size;
    }

    /// \brief Difference  operator
    bool operator!=(const OrthantConeTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Returns the dimension of the box.
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

    /// \brief Check whether a vector x lies within the orthant.
    ///
    /// \param[in] x vector to check .
    ///
    template<typename VectorLike>
    bool isInside(const Eigen::MatrixBase<VectorLike> & x, const Scalar prec = Scalar(0)) const
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
    void project(
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

    using Base::derived;
    /// \brief Returns the dual cone associated with this.
    ///
    /// \remarks Orthant cone are by definition self dual.
    DualCone dual() const
    {
      return derived();
    }

  protected:
    Eigen::DenseIndex m_size;
  }; // struct OrthantConeTpl

} // namespace pinocchio
#endif // ifndef __pinocchio_algorithm_constraints_orthant_cone_hpp__
