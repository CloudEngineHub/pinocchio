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
  const Eigen::DenseIndex rows = 10, cols = 20, inner_dim = 30;
  for (int i = 0; i < N; ++i)
  {
    const Eigen::MatrixXd lhs = Eigen::MatrixXd::Random(rows, inner_dim),
                          rhs = Eigen::MatrixXd::Random(inner_dim, cols);
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(rows, cols);
    const auto res_gt = (lhs * rhs).eval();
    matrix_product(lhs, rhs, res);

    BOOST_CHECK(res.isApprox(res_gt, 1e-14));
  }
}

BOOST_AUTO_TEST_SUITE_END()
