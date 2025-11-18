//
// Copyright (c) 2025 INRIA
//

#include <iostream>

#include "pinocchio/utils/eigen.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

template<typename Base>
struct Accessor : Base
{
  using Base::is_partial_static_size_product;
  using Base::is_static_size_product;
};

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_helpers)
{
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 2);

  BOOST_CHECK(!helper::is_eigen_noalias_v<decltype(A)>);
  BOOST_CHECK(helper::is_eigen_noalias_v<decltype(A.noalias())>);

  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix3d;
  BOOST_CHECK(helper::has_fixed_rows_v<Matrix3d>);
  BOOST_CHECK(!helper::has_fixed_cols_v<Matrix3d>);
  BOOST_CHECK(!helper::has_fixed_size_v<Matrix3d>);

  typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Matrixd3;
  BOOST_CHECK(!helper::has_fixed_rows_v<Matrixd3>);
  BOOST_CHECK(helper::has_fixed_cols_v<Matrixd3>);
  BOOST_CHECK(!helper::has_fixed_size_v<Matrixd3>);

  typedef Eigen::Matrix<double, 3, 3> Matrix33;
  BOOST_CHECK(helper::has_fixed_rows_v<Matrix33>);
  BOOST_CHECK(helper::has_fixed_cols_v<Matrix33>);
  BOOST_CHECK(helper::has_fixed_size_v<Matrix33>);
}

BOOST_AUTO_TEST_CASE(test_dynamic_matrix)
{
  const Eigen::DenseIndex n = 10, m = 5;
  Eigen::MatrixXd A = Eigen::MatrixXd::Constant(n, n, 1);
  Eigen::MatrixXd B = Eigen::MatrixXd::Constant(n, m, 2);
  Eigen::MatrixXd C = Eigen::MatrixXd::Random(n, m);

  const auto C_expression = A * B;

  BOOST_CHECK(C_expression.rows() == A.rows());
  BOOST_CHECK(C_expression.cols() == B.cols());

  auto C_op = promote_static_op(C);
  typedef decltype(C_op) PromotedType;
  BOOST_CHECK(&C_op.expression() == &C);

  BOOST_CHECK(!Accessor<PromotedType>::is_static_size_product<decltype(C_expression)>());
  BOOST_CHECK(!Accessor<PromotedType>::is_partial_static_size_product<decltype(C_expression)>());

  C_op = A * B;
  const auto res_aliasing = C_expression.eval();
  BOOST_CHECK(C == res_aliasing);

  // Test with noalias
  A.setConstant(3);
  B.setConstant(4);

  auto C_noalias_op = promote_static_op(C.noalias());
  typedef decltype(C_noalias_op) NoAliasPromotedType;
  BOOST_CHECK(&C_noalias_op.expression().expression() == &C);

  BOOST_CHECK(!Accessor<NoAliasPromotedType>::is_static_size_product<decltype(C_expression)>());
  BOOST_CHECK(
    !Accessor<NoAliasPromotedType>::is_partial_static_size_product<decltype(C_expression)>());

  C_noalias_op = A * B;
  const auto res_noaliasing = C_expression.eval();
  BOOST_CHECK(res_noaliasing != res_aliasing);
  BOOST_CHECK(C == res_noaliasing);
}

BOOST_AUTO_TEST_CASE(test_static_matrix)
{
  constexpr Eigen::DenseIndex Rows = 10, Cols = 5, InnerDim = 8;
  typedef Eigen::Matrix<double, Rows, InnerDim> LhsType;
  typedef Eigen::Matrix<double, InnerDim, Cols> RhsType;
  typedef Eigen::Matrix<double, Rows, Cols> ResType;

  LhsType A = LhsType::Constant(1);
  RhsType B = RhsType::Constant(2);
  ResType C = ResType::Random();

  const auto C_expression = A * B;

  BOOST_CHECK(C_expression.rows() == A.rows());
  BOOST_CHECK(C_expression.cols() == B.cols());

  auto C_op = promote_static_op(C);
  typedef decltype(C_op) PromotedType;
  BOOST_CHECK(&C_op.expression() == &C);

  BOOST_CHECK(Accessor<PromotedType>::is_static_size_product<decltype(C_expression)>());
  BOOST_CHECK(Accessor<PromotedType>::is_partial_static_size_product<decltype(C_expression)>());

  C_op = A * B;
  const auto res_aliasing = C_expression.eval();
  BOOST_CHECK(C == res_aliasing);

  // Test with noalias
  A.setConstant(3);
  B.setConstant(4);

  auto C_noalias_op = promote_static_op(C.noalias());
  typedef decltype(C_noalias_op) NoAliasPromotedType;
  BOOST_CHECK(&C_noalias_op.expression().expression() == &C);

  BOOST_CHECK(Accessor<NoAliasPromotedType>::is_static_size_product<decltype(C_expression)>());
  BOOST_CHECK(
    Accessor<NoAliasPromotedType>::is_partial_static_size_product<decltype(C_expression)>());

  C_noalias_op = A * B;
  const auto res_noaliasing = C_expression.eval();
  BOOST_CHECK(res_noaliasing != res_aliasing);
  BOOST_CHECK(C == res_noaliasing);
}

BOOST_AUTO_TEST_SUITE_END()
