//
// Copyright (c) 2025 INRIA, KU Leuven
//

#ifndef __pinocchio_algorithm_ip_solver_cone_operations_hpp__
#define __pinocchio_algorithm_ip_solver_cone_operations_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/utils/reference.hpp"

namespace pinocchio
{
  template<typename _Scalar>
  struct IPMSolverConeOperations
  {
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef Eigen::Matrix<Scalar, 2, 1> Vector2;
    typedef Eigen::Matrix<Scalar, 2, 2> Matrix2x2;
    typedef Eigen::Matrix<Scalar, 3, 3> Matrix3x3;
    typedef CoulombFrictionConeTpl<Scalar> Cone;

    /// \brief compute the cone product of two 3-vectors x and y
    template<typename Vector3Like>
    static Vector3
    coneProduct(const Eigen::MatrixBase<Vector3Like> & x, const Eigen::MatrixBase<Vector3Like> & y);

    template<typename Vector3Like>
    static void coneProduct(
      const Eigen::MatrixBase<Vector3Like> & x,
      const Eigen::MatrixBase<Vector3Like> & y,
      Eigen::MatrixBase<Vector3Like> & out);

    /// \brief compute the inverse cone product of two 3-vectors x and y
    template<typename Vector3Like>
    static Vector3 inverseConeProduct(
      const Eigen::MatrixBase<Vector3Like> & x, const Eigen::MatrixBase<Vector3Like> & y);

    /// \brief compute the inverse cone product of two concatenated lists 3-vectors x and y, store
    /// the result in out
    template<typename Vector3Like>
    static void inverseConeProduct(
      const Eigen::MatrixBase<Vector3Like> & x,
      const Eigen::MatrixBase<Vector3Like> & y,
      Eigen::MatrixBase<Vector3Like> & out);

    /// \brief evaluates x = H(\lambda^0.5) * x
    template<typename Vector3Like>
    static Vector3
    scale2(const Eigen::MatrixBase<Vector3Like> & lambda, const Eigen::MatrixBase<Vector3Like> & x);

    /// \brief evaluates x = H(\lambda^{-0.5}) * x
    template<typename Vector3Like>
    static Vector3 scale2Inv(
      const Eigen::MatrixBase<Vector3Like> & lambda, const Eigen::MatrixBase<Vector3Like> & x);

    /// \brief evaluates out = H(\lambda^0.5) * x, for concatenated lists of 3-vectors
    template<typename VectorLike>
    static void scale2(
      const Eigen::MatrixBase<VectorLike> & lambda,
      const Eigen::MatrixBase<VectorLike> & x,
      Eigen::MatrixBase<VectorLike> & out);

    /// \brief evaluates out = H(\lambda^{-0.5}) * x, for concatenated lists of 3-vectors
    template<typename VectorLike>
    static void scale2Inv(
      const Eigen::MatrixBase<VectorLike> & lambda,
      const Eigen::MatrixBase<VectorLike> & x,
      Eigen::MatrixBase<VectorLike> & out);

    /// \brief represents a scaling matrix of a pair (s, z) of 3-vectors
    struct ScalingMatrix
    {
      /// \brief computes the scaling matrix for the pair (s, z)
      Vector3 compute(const Vector3 & s, const Vector3 & z);

      /// \brief applies the scaling matrix to a 3-vector x, res = W*x
      Vector3 apply(const Vector3 & x);

      /// \brief applies the inverse scaling matrix to a 3-vector x, res = W^{-1}*x
      Vector3 applyInverse(const Vector3 & x);

      /// \brief returns the scaling matrix
      Matrix3x3 getMatrix();

      /// \brief returns the inverse scaling matrix
      Matrix3x3 getInverseMatrix() const;

      Vector3 v;
      Scalar beta;
    };

  protected:
    /// \brief compute x^T J y
    /// Note:  J = [-I 0]
    ///            [ 0 1]
    template<typename Vector3Like>
    static Scalar
    jdot(const Eigen::MatrixBase<Vector3Like> & x, const Eigen::MatrixBase<Vector3Like> & y);

    /// \brief compute Jx
    template<typename Vector3Like>
    static Vector3 jprod(const Eigen::MatrixBase<Vector3Like> & x);

    /// \brief Compute the J norm: sqrt(x^T J x)
    template<typename Vector3Like>
    static Scalar jnorm(const Eigen::MatrixBase<Vector3Like> & x);

    /// \brief Compute the J squared norm: x^T J x
    template<typename Vector3Like>
    static Scalar jnorm2(const Eigen::MatrixBase<Vector3Like> & x);
  };

  // implementation

  template<typename Scalar>
  template<typename Vector3Like>
  Scalar IPMSolverConeOperations<Scalar>::jdot(
    const Eigen::MatrixBase<Vector3Like> & x, const Eigen::MatrixBase<Vector3Like> & y)
  {
    Scalar xn = x.coeff(2);
    Scalar yn = y.coeff(2);
    auto xt = x.template head<2>();
    auto yt = y.template head<2>();

    return xn * yn - xt.dot(yt);
  }

  template<typename Scalar>
  template<typename Vector3Like>
  Eigen::Matrix<Scalar, 3, 1>
  IPMSolverConeOperations<Scalar>::jprod(const Eigen::MatrixBase<Vector3Like> & x)
  {
    return Vector3(-x[0], -x[1], x[2]);
  }

  template<typename Scalar>
  template<typename Vector3Like>
  Eigen::Matrix<Scalar, 3, 1> IPMSolverConeOperations<Scalar>::scale2(
    const Eigen::MatrixBase<Vector3Like> & lambda, const Eigen::MatrixBase<Vector3Like> & x)
  {
    Vector3 ret;
    Scalar a = jnorm(lambda);
    Scalar lx = jdot(lambda, x) / a;
    Scalar xN = x[2];
    Scalar c = -(lx + xN) / (lambda[2] / a + 1) / a;
    ret.template head<2>() = x.template head<2>() + c * lambda.template head<2>();
    ret[2] = lx;
    ret *= 1.0 / a;
    return ret;
  }

  template<typename Scalar>
  template<typename Vector3Like>
  Eigen::Matrix<Scalar, 3, 1> IPMSolverConeOperations<Scalar>::scale2Inv(
    const Eigen::MatrixBase<Vector3Like> & lambda, const Eigen::MatrixBase<Vector3Like> & x)
  {
    Vector3 ret;
    Scalar a = jnorm(lambda);
    Scalar lx = lambda.dot(x) / a;
    Scalar xN = x[2];
    Scalar c = (lx + xN) / (lambda[2] / a + 1) / a;
    ret.template head<2>() = x.template head<2>() + c * lambda.template head<2>();
    ret[2] = lx;
    ret *= a;
    return ret;
  }

  template<typename Scalar>
  Eigen::Matrix<Scalar, 3, 1>
  IPMSolverConeOperations<Scalar>::ScalingMatrix::compute(const Vector3 & s, const Vector3 & z)
  {
    Vector3 lambda;
    Scalar aa = jnorm(s);
    Scalar bb = jnorm(z);
    beta = math::sqrt(aa / bb);
    Vector3 s_stripe = s / aa;
    Vector3 z_stripe = z / bb;
    double cc = math::sqrt((1.0 + s.dot(z) / aa / bb) / 2.);
    v = 1. / (2 * cc) * (s_stripe + jprod(z_stripe));
    v[2] += 1.;
    v *= 1 / math::sqrt(2.0 * v[2]);
    lambda[2] = cc;
    Scalar dd = 2 * cc + s_stripe[2] + z_stripe[2];
    lambda.template head<2>() = (cc + z_stripe[2]) / dd * s_stripe.template head<2>()
                                + (cc + s_stripe[2]) / dd * z_stripe.template head<2>();
    lambda *= math::sqrt(aa * bb);
    return lambda;
  }

  template<typename Scalar>
  Eigen::Matrix<Scalar, 3, 1>
  IPMSolverConeOperations<Scalar>::ScalingMatrix::apply(const Vector3 & x)
  {
    Scalar w = x.dot(v);
    return beta * (2 * v * w - jprod(x));
  };

  template<typename Scalar>
  Eigen::Matrix<Scalar, 3, 1>
  IPMSolverConeOperations<Scalar>::ScalingMatrix::applyInverse(const Vector3 & x)
  {
    return -1. / beta * (jprod(2 * v * (jprod(-x).dot(v)) + x));
  };

  template<typename Scalar>
  Eigen::Matrix<Scalar, 3, 3> IPMSolverConeOperations<Scalar>::ScalingMatrix::getMatrix()
  {
    Matrix3x3 ret = Matrix3x3::Zero();
    ret = 2 * v * v.transpose();
    ret.template block<2, 2>(0, 0) += Matrix2x2::Identity();
    ret(2, 2) -= 1;
    return beta * ret;
  };

  template<typename Scalar>
  Eigen::Matrix<Scalar, 3, 3>
  IPMSolverConeOperations<Scalar>::ScalingMatrix::getInverseMatrix() const
  {
    Matrix3x3 ret = Matrix3x3::Zero();
    ret = 2 * v * (-jprod(v)).transpose();
    ret += Matrix3x3::Identity();
    ret.row(2) *= -1;
    return 1.0 / beta * ret;
  };

  template<typename Scalar>
  template<typename VectorLike>
  void IPMSolverConeOperations<Scalar>::scale2(
    const Eigen::MatrixBase<VectorLike> & lambda,
    const Eigen::MatrixBase<VectorLike> & x,
    Eigen::MatrixBase<VectorLike> & out)
  {
    for (int i = 0; i < x.size(); i += 3)
    {
      out.template segment<3>(i) = scale2(lambda.template segment<3>(i), x.template segment<3>(i));
    }
  }

  template<typename Scalar>
  template<typename VectorLike>
  void IPMSolverConeOperations<Scalar>::scale2Inv(
    const Eigen::MatrixBase<VectorLike> & lambda,
    const Eigen::MatrixBase<VectorLike> & x,
    Eigen::MatrixBase<VectorLike> & out)
  {
    for (int i = 0; i < x.size(); i += 3)
    {
      out.template segment<3>(i) =
        scale2Inv(lambda.template segment<3>(i), x.template segment<3>(i));
    }
  }

  template<typename Scalar>
  template<typename Vector3Like>
  Eigen::Matrix<Scalar, 3, 1> IPMSolverConeOperations<Scalar>::coneProduct(
    const Eigen::MatrixBase<Vector3Like> & x, const Eigen::MatrixBase<Vector3Like> & y)
  {
    Vector3 ret;
    ret.template head<2>() = y[2] * x.template head<2>() + x[2] * y.template head<2>();
    ret[2] = x.dot(y);
    return ret;
  }

  template<typename Scalar>
  template<typename Vector3Like>
  void IPMSolverConeOperations<Scalar>::coneProduct(
    const Eigen::MatrixBase<Vector3Like> & x,
    const Eigen::MatrixBase<Vector3Like> & y,
    Eigen::MatrixBase<Vector3Like> & out)
  {
    for (int i = 0; i < x.size(); i += 3)
    {
      out.template segment<3>(i) = coneProduct(x.template segment<3>(i), y.template segment<3>(i));
    }
  };

  template<typename Scalar>
  template<typename Vector3Like>
  Eigen::Matrix<Scalar, 3, 1> IPMSolverConeOperations<Scalar>::inverseConeProduct(
    const Eigen::MatrixBase<Vector3Like> & x, const Eigen::MatrixBase<Vector3Like> & y)
  {
    // see ecos implementation paper
    Vector3 ret;
    Scalar q = jnorm2(x);
    Scalar nu = x.template head<2>().dot(y.template head<2>());
    ret.template head<2>() =
      1. / q * (nu / x[2] - y[2]) * x.template head<2>() + 1. / x[2] * y.template head<2>();
    ret[2] = 1. / q * (x[2] * y[2] - nu);
    return ret;
  }

  template<typename Scalar>
  template<typename Vector3Like>
  void IPMSolverConeOperations<Scalar>::inverseConeProduct(
    const Eigen::MatrixBase<Vector3Like> & x,
    const Eigen::MatrixBase<Vector3Like> & y,
    Eigen::MatrixBase<Vector3Like> & out)
  {
    for (Eigen::Index i = 0; i < x.size(); i += 3)
    {
      out.template segment<3>(i) =
        inverseConeProduct(x.template segment<3>(i), y.template segment<3>(i));
    }
  };

  template<typename Scalar>
  template<typename Vector3Like>
  Scalar IPMSolverConeOperations<Scalar>::jnorm(const Eigen::MatrixBase<Vector3Like> & x)
  {
    return math::sqrt(jnorm2(x));
  }

  template<typename Scalar>
  template<typename Vector3Like>
  Scalar IPMSolverConeOperations<Scalar>::jnorm2(const Eigen::MatrixBase<Vector3Like> & x)
  {
    Scalar xn = x.coeff(2);
    auto xt = x.template head<2>();

    Scalar xn2 = xn * xn;
    Scalar xt2 = xt.dot(xt);

    assert(xn2 - xt2 >= 0);
    return xn2 - xt2;
  }

} // namespace pinocchio

#endif // __pinocchio_algorithm_ip_solver_cone_operations_hpp__
