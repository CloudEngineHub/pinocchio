//
// Copyright (c) 2026 INRIA
//

#include <limits>
#include <pinocchio/math.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

using namespace pinocchio;
typedef Eigen::Matrix<double, 1, 1> M11;
typedef BlockDiagonalMatrix::Matrix Matrix;
typedef BlockDiagonalMatrix::MatrixMap MatrixMap;
typedef BlockDiagonalMatrix::ConstMatrixMap ConstMatrixMap;
typedef BlockDiagonalMatrix::Vector Vector;
typedef MatrixBlockElementTpl<Matrix> MatrixBlockElement;
typedef MatrixBlockElementTpl<MatrixMap> MatrixMapBlockElement;
typedef MatrixBlockElementTpl<const Matrix> ConstMatrixBlockElement;
typedef MatrixBlockElementTpl<ConstMatrixMap> ConstMatrixMapBlockElement;

template<typename _MatrixBlockElement>
void test_assignment(const MatrixBlockElementPlain<_MatrixBlockElement> & matrix_block_element)
{
  const auto size = matrix_block_element.size();
  const Matrix square_matrix = Matrix::Random(size, size);

  const auto matrix_block_element_plain = matrix_block_element.matrix();

  const int num_tests =
#ifdef NDEBUG
    100000
#else
    1000
#endif
    ;

  // evalTo
  for (int k = 0; k < num_tests; ++k)
  {
    Matrix res = Matrix::Random(size, size);
    matrix_block_element.evalTo(res);

    BOOST_CHECK(res.isApprox(matrix_block_element_plain));
  }

  // addTo
  for (int k = 0; k < num_tests; ++k)
  {
    const Matrix res_ref = square_matrix + matrix_block_element_plain;
    Matrix res(square_matrix);
    matrix_block_element.addTo(res);

    BOOST_CHECK(res.isApprox(res_ref));
  }

  // subTo
  for (int k = 0; k < num_tests; ++k)
  {
    const Matrix res_ref = square_matrix - matrix_block_element_plain;
    Matrix res(square_matrix);
    matrix_block_element.subTo(res);

    BOOST_CHECK(res.isApprox(res_ref));
  }
}

template<typename _MatrixBlockElement>
void test_inverse(const MatrixBlockElementPlain<_MatrixBlockElement> & matrix_block_element)
{
  const auto size = matrix_block_element.size();
  const auto matrix_block_element_plain = matrix_block_element.matrix();
  const auto matrix_block_element_inverse = matrix_block_element.inverse();
  const auto matrix_block_element_inverse_plain = matrix_block_element_inverse.matrix();

  if (matrix_block_element.type() == pinocchio::MatrixBlockType::Zero)
  {
    BOOST_CHECK(
      matrix_block_element_inverse_plain
      == Eigen::MatrixXd::Constant(size, size, std::numeric_limits<double>::infinity()));
  }
  else
  {
    BOOST_CHECK(matrix_block_element_plain.inverse().isApprox(matrix_block_element_inverse_plain));
  }

  auto matrix_block_element_copy = matrix_block_element.derived();

  const int num_tests =
#ifdef NDEBUG
    100000
#else
    1000
#endif
    ;

  for (int k = 0; k < num_tests; ++k)
  {
    matrix_block_element_copy.setRandom();
    const auto matrix_block_element_copy_inverse = matrix_block_element_copy.inverse();

    if (matrix_block_element_copy.type() == pinocchio::MatrixBlockType::Zero)
    {
      BOOST_CHECK(
        matrix_block_element_copy_inverse.matrix()
        == Eigen::MatrixXd::Constant(size, size, std::numeric_limits<double>::infinity()));
    }
    else
    {
      BOOST_CHECK(matrix_block_element_copy.matrix().inverse().isApprox(
        matrix_block_element_copy_inverse.matrix()));
    }
  }
}

template<typename _MatrixBlockElement>
void test_operations(const MatrixBlockElementPlain<_MatrixBlockElement> & _matrix_block_element)
{
  const auto & matrix_block_element = _matrix_block_element.derived();

  const auto size = matrix_block_element.size();
  const Vector diagonal_vector = Vector::Random(size);

  const auto matrix_block_element_plain = matrix_block_element.matrix();

  const int num_tests =
#ifdef NDEBUG
    100000
#else
    1000
#endif
    ;

  // operator+
  for (int k = 0; k < num_tests; ++k)
  {
    auto block_sum = matrix_block_element + diagonal_vector.asDiagonal();

    MatrixBlockElement res_block = {pinocchio::MatrixBlockType::Diagonal, size, diagonal_vector};
    res_block = block_sum;

    const auto res_eigen =
      (matrix_block_element.matrix() + Matrix(diagonal_vector.asDiagonal())).eval();

    BOOST_CHECK(Matrix(res_block.container().asDiagonal()).isApprox(res_eigen));
  }

  // operator-
  for (int k = 0; k < num_tests; ++k)
  {
    auto block_sum = matrix_block_element - diagonal_vector.asDiagonal();

    MatrixBlockElement res_block = {pinocchio::MatrixBlockType::Diagonal, size, diagonal_vector};
    res_block = block_sum;

    const auto res_eigen =
      (matrix_block_element.matrix() - Matrix(diagonal_vector.asDiagonal())).eval();

    BOOST_CHECK(Matrix(res_block.container().asDiagonal()).isApprox(res_eigen));
  }
}

// void test_applyOnTheRight(const BlockDiagonalMatrix & block_diagonal_matrix)
// {
//   const auto rows = block_diagonal_matrix.rows();
//   const auto cols = 20;

//   const auto bdm_plain = block_diagonal_matrix.matrix();

//   const int num_tests =
// #ifdef NDEBUG
//     100000
// #else
//     1000
// #endif
//     ;

//   // assign
//   for (int k = 0; k < num_tests; ++k)
//   {
//     const Matrix rhs_matrix = Matrix::Random(rows, cols);
//     const auto res_ref = (bdm_plain * rhs_matrix).eval();
//     const auto res = block_diagonal_matrix * rhs_matrix;
//     BOOST_CHECK(res.isApprox(res_ref));
//   }

//   // add_assign
//   for (int k = 0; k < num_tests; ++k)
//   {
//     const Matrix rhs_matrix = Matrix::Random(rows, cols);
//     const Matrix res0 = Matrix::Random(rows, cols);
//     const auto res_ref = (res0 + bdm_plain * rhs_matrix).eval();
//     Matrix res = res0;
//     block_diagonal_matrix.applyOnTheRight<pinocchio::internal::add_assign_op>(rhs_matrix, res);
//     BOOST_CHECK(res.isApprox(res_ref));
//   }

//   // sub_assign
//   for (int k = 0; k < num_tests; ++k)
//   {
//     const Matrix rhs_matrix = Matrix::Random(rows, cols);
//     const Matrix res0 = Matrix::Random(rows, cols);
//     const auto res_ref = (res0 - bdm_plain * rhs_matrix).eval();
//     Matrix res = res0;
//     block_diagonal_matrix.applyOnTheRight<pinocchio::internal::sub_assign_op>(rhs_matrix, res);
//     BOOST_CHECK(res.isApprox(res_ref));
//   }
// }

BOOST_AUTO_TEST_CASE(test_matrix_block_elements_eigen_map)
{
  typedef MatrixBlockElementTpl<MatrixMap> MatrixBlockElement;
  typedef MatrixBlockElementTpl<ConstMatrixMap> ConstMatrixBlockElement;
  const Eigen::Index size = 10;

  // Test non const version
  {
    Matrix diagonal_vector = Matrix::Ones(size, 1);

    const auto matrix_map = make_map<MatrixMap>(diagonal_vector);
    MatrixBlockElement lhs_block = {pinocchio::MatrixBlockType::Diagonal, size, matrix_map};

    BOOST_CHECK(lhs_block.map == diagonal_vector);
    BOOST_CHECK(lhs_block.map.data() == diagonal_vector.data());

    test_assignment(lhs_block);
    test_operations(lhs_block);
  }

  // Test const version
  {
    Matrix diagonal_vector = Matrix::Ones(size, 1);

    const auto matrix_map = make_map<ConstMatrixMap>(diagonal_vector);
    ConstMatrixBlockElement lhs_block = {pinocchio::MatrixBlockType::Diagonal, size, matrix_map};

    BOOST_CHECK(lhs_block.map == diagonal_vector);
    BOOST_CHECK(lhs_block.map.data() == diagonal_vector.data());

    test_assignment(lhs_block);
  }
}

BOOST_AUTO_TEST_CASE(test_matrix_block_elements_eigen_matrix)
{
  const Eigen::Index size = 10;
  typedef MatrixBlockElementTpl<Matrix> MatrixBlockElement;
  typedef MatrixBlockElementTpl<const Matrix> ConstMatrixBlockElement;

  // Test non const version
  {
    Matrix diagonal_vector = Matrix::Ones(size, 1);

    MatrixBlockElement lhs_block = {pinocchio::MatrixBlockType::Diagonal, size, diagonal_vector};

    BOOST_CHECK(lhs_block.container() == diagonal_vector);

    test_assignment(lhs_block);
    test_operations(lhs_block);
  }

  // Test const version
  {
    const Matrix diagonal_vector = Matrix::Ones(size, 1);

    ConstMatrixBlockElement lhs_block = {
      pinocchio::MatrixBlockType::Diagonal, size, diagonal_vector};

    BOOST_CHECK(lhs_block.container() == diagonal_vector);

    test_assignment(lhs_block);
  }
}

BOOST_AUTO_TEST_CASE(test_inverse_method)
{
  const Eigen::Index size = 10;
  typedef MatrixBlockElementTpl<Matrix> MatrixBlockElement;

  // Zero block
  {
    MatrixBlockElement matrix_block_element = {pinocchio::MatrixBlockType::Zero, size};

    test_inverse(matrix_block_element);
  }

  // Identity block
  {
    MatrixBlockElement matrix_block_element = {pinocchio::MatrixBlockType::Identity, size};

    test_inverse(matrix_block_element);
  }

  // Scalar Identity block
  {
    const double scale_value = 1.;
    M11 scale_mat = M11(scale_value);
    MatrixBlockElement matrix_block_element = {
      pinocchio::MatrixBlockType::ScalarIdentity, size, scale_mat};

    test_inverse(matrix_block_element);
  }

  // Diagonal
  {
    Matrix diagonal_vector = Matrix::Ones(size, 1);
    MatrixBlockElement matrix_block_element = {
      pinocchio::MatrixBlockType::Diagonal, size, diagonal_vector};

    test_inverse(matrix_block_element);
  }

  // Plain
  {
    Matrix plain_matrix = Matrix::Identity(size, size);
    MatrixBlockElement matrix_block_element = {
      pinocchio::MatrixBlockType::Plain, size, plain_matrix};

    test_inverse(matrix_block_element);
  }
}

BOOST_AUTO_TEST_SUITE_END()
