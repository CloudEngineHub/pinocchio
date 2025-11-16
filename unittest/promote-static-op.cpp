//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/utils/eigen.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

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
  BOOST_CHECK(&C_op.expression() == &C);

  C_op = A * B;
  const auto res_aliasing = C_expression.eval();
  BOOST_CHECK(C == res_aliasing);

  // Test with noalias
  A.setConstant(3);
  B.setConstant(4);

  auto C_noalias_op = promote_static_op(C.noalias());
  BOOST_CHECK(&C_noalias_op.expression().expression() == &C);

  C_noalias_op = A * B;
  const auto res_noaliasing = C_expression.eval();
  BOOST_CHECK(res_noaliasing != res_aliasing);
  BOOST_CHECK(C == res_noaliasing);
}

BOOST_AUTO_TEST_SUITE_END()
