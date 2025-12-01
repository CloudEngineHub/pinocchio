//
// Copyright (c) 2025 INRIA
//

#include <iostream>

#include <pinocchio/math/matrix-product.hpp>

#include <boost/variant.hpp> // to avoid C99 warnings

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;

#ifndef NDEBUG
const int N = int(1e3);
#else
const int N = int(1e4);
#endif

BOOST_AUTO_TEST_CASE(test_col_major)
{
  typedef Eigen::MatrixXd Matrix;
  const Eigen::DenseIndex rows = 10, cols = 20, inner_dim = 30;
  for (int i = 0; i < N; ++i)
  {
    const Matrix lhs = Matrix::Random(rows, inner_dim), rhs = Matrix::Random(inner_dim, cols);
    Matrix res = Matrix::Zero(rows, cols);
    const auto res_gt = (lhs * rhs).eval();
    matrix_product(lhs, rhs, res);

    BOOST_CHECK(res.isApprox(res_gt, 1e-14));
  }
}

BOOST_AUTO_TEST_CASE(test_row_major)
{
  typedef PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(Eigen::MatrixXd) Matrix;

  const Eigen::DenseIndex rows = 10, cols = 20, inner_dim = 30;
  for (int i = 0; i < N; ++i)
  {
    const Matrix lhs = Matrix::Random(rows, inner_dim), rhs = Matrix::Random(inner_dim, cols);
    Matrix res = Matrix::Zero(rows, cols);
    const auto res_gt = (lhs * rhs).eval();
    matrix_product(lhs, rhs, res);

    BOOST_CHECK(res.isApprox(res_gt, 1e-14));
  }
}

BOOST_AUTO_TEST_SUITE_END()
