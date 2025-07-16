//
// Copyright (c) 2025 INRIA
//

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

#include "pinocchio/container/matrix-stack.hpp"
// #define ALIGNMENT_VALUE 128
#define ALIGNMENT_VALUE 8

using namespace pinocchio;
typedef Eigen::MatrixXf MatrixXs;
typedef MatrixXs::Scalar Scalar;
typedef MatrixStackTpl<MatrixXs, ALIGNMENT_VALUE> MatrixXsStack;
typedef MatrixStackTpl<PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXs), ALIGNMENT_VALUE>
  RowMatrixXsStack;
// typedef EigenStorageTpl<Eigen::VectorXd> EigenStorageVector;

bool is_aligned(const void * ptr, const std::size_t alignment)
{
  assert(
    alignment >= sizeof(void *) && (alignment & (alignment - 1)) == 0
    && "Alignment must be at least sizeof(void*) and a power of 2");
  return reinterpret_cast<std::size_t>(ptr) % alignment == 0;
}

template<typename T>
const T & as_const(T & value)
{
  return value;
}

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(print_info)
{
  std::cout << "sizeof(Scalar): " << sizeof(Scalar) << std::endl;
  std::cout << "alignment: " << MatrixXsStack::Alignment << std::endl;
  std::cout << "sizeof(void*): " << sizeof(void *) << std::endl;
}

BOOST_AUTO_TEST_CASE(matrix_stack_empty)
{
  MatrixXsStack matrix_stack(0);
  BOOST_CHECK(matrix_stack.size() == 0);
  BOOST_CHECK(matrix_stack.capacity() == 0);
  BOOST_CHECK(matrix_stack.data() == nullptr);

  MatrixXsStack matrix_stack_copy(matrix_stack);
  BOOST_CHECK(matrix_stack_copy.size() == 0);
  BOOST_CHECK(matrix_stack_copy.capacity() == 0);
  BOOST_CHECK(matrix_stack_copy.data() == nullptr);
}

BOOST_AUTO_TEST_CASE(matrix_stack_default)
{

  BOOST_CHECK(ALIGNMENT_VALUE == MatrixXsStack::Alignment);

  { // case with no allocation
    MatrixXsStack matrix_stack(100);
    BOOST_CHECK(matrix_stack.size() == 0);
    BOOST_CHECK(matrix_stack.capacity() == 100);

    BOOST_CHECK(matrix_stack.data() == nullptr);
    BOOST_CHECK(is_aligned(matrix_stack.data(), MatrixXsStack::Alignment));

    // First push_back
    matrix_stack.push_back(2, 2);
    BOOST_CHECK(matrix_stack.data() != nullptr);
    BOOST_CHECK(matrix_stack.size() == 1);
    BOOST_CHECK(is_aligned(matrix_stack.data(), MatrixXsStack::Alignment));
    BOOST_CHECK(matrix_stack.back().rows() == 2);
    BOOST_CHECK(matrix_stack.back().cols() == 2);
    BOOST_CHECK(matrix_stack.data() == matrix_stack.back().data());

    matrix_stack.back().setOnes();
    BOOST_CHECK(matrix_stack.back().isOnes(0));

    // Second push_back
    matrix_stack.push_back(3, 3);
    BOOST_CHECK(matrix_stack.size() == 2);
    BOOST_CHECK(matrix_stack.data() != matrix_stack.back().data());
    BOOST_CHECK(is_aligned(matrix_stack.back().data(), MatrixXsStack::Alignment));
    BOOST_CHECK(matrix_stack.back().rows() == 3);
    BOOST_CHECK(matrix_stack.back().cols() == 3);
    matrix_stack.back().setConstant(2);

    BOOST_CHECK(matrix_stack[0].isOnes(0));
    BOOST_CHECK(matrix_stack[1].isConstant(2, 0));

    // operator==
    BOOST_CHECK(matrix_stack == matrix_stack);

    // Test copy constructor
    MatrixXsStack matrix_stack_copy = matrix_stack;
    BOOST_CHECK(matrix_stack_copy == matrix_stack);

    matrix_stack_copy.push_back(1, 6);
    BOOST_CHECK(matrix_stack_copy != matrix_stack);
  }

  { // case with memory allocation
    const Eigen::DenseIndex max_rows = 3;
    const Eigen::DenseIndex max_cols = 4;
    const std::size_t matrix_max_size = max_rows * max_cols;
    MatrixXsStack matrix_stack(100, matrix_max_size);
    BOOST_CHECK(matrix_stack.size() == 0);
    BOOST_CHECK(matrix_stack.capacity() == 100);

    BOOST_CHECK(matrix_stack.data() != 0);
    BOOST_CHECK(is_aligned(matrix_stack.data(), MatrixXsStack::Alignment));

    // First push_back
    matrix_stack.push_back(2, 2);
    BOOST_CHECK(matrix_stack.data() == matrix_stack.back().data());
    BOOST_CHECK(is_aligned(matrix_stack.back().data(), MatrixXsStack::Alignment));
    BOOST_CHECK(matrix_stack.back().rows() == 2);
    BOOST_CHECK(matrix_stack.back().cols() == 2);
    matrix_stack.back().setOnes();

    // Check constness
    BOOST_CHECK(as_const(matrix_stack).data() == as_const(matrix_stack).back().data());
    BOOST_CHECK(as_const(matrix_stack).back().rows() == 2);
    BOOST_CHECK(as_const(matrix_stack).back().cols() == 2);
    BOOST_CHECK(as_const(matrix_stack).back().isOnes(0.));

    // Check a second push_back
    matrix_stack.push_back(3, 3);
    BOOST_CHECK(matrix_stack.size() == 2);
    BOOST_CHECK(matrix_stack.data() != matrix_stack.back().data());
    BOOST_CHECK(is_aligned(matrix_stack.back().data(), MatrixXsStack::Alignment));
    BOOST_CHECK(matrix_stack.back().rows() == 3);
    BOOST_CHECK(matrix_stack.back().cols() == 3);
    matrix_stack.back().setConstant(2);

    // Check element values
    BOOST_CHECK(matrix_stack[0].isOnes(0.));
    BOOST_CHECK(matrix_stack[1].isConstant(2., 0.));

    // Check distance between the two pointers
    const std::size_t first_elt_raw_size =
      std::size_t(matrix_stack[0].size()) * sizeof(MatrixXsStack::Scalar);
    const std::size_t max_distance = first_elt_raw_size + MatrixXsStack::Alignment;
    const std::size_t ptr_distance = reinterpret_cast<std::size_t>(matrix_stack[1].data())
                                     - reinterpret_cast<std::size_t>(matrix_stack[0].data());

    BOOST_CHECK(ptr_distance <= max_distance);
    BOOST_CHECK(ptr_distance >= first_elt_raw_size);

    std::cout << "ptr_distance: " << ptr_distance << std::endl;
    std::cout << "first_elt_raw_size: " << first_elt_raw_size << std::endl;
    std::cout << "max_distance: " << max_distance << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
