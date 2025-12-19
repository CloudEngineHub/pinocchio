//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/utils/size-in-bytes.hpp"
#include "pinocchio/utils/std-array.hpp"
#include "pinocchio/utils/std-vector.hpp"

#include <cstddef>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

template<int size>
struct SimpleStruct
{
  SimpleStruct() = default;

  std::size_t sizeInBytes() const
  {
    return sizeof(char) * size;
  }

  char data[size];
};

typedef SimpleStruct<1> SimpleStruct1;
typedef SimpleStruct<10> SimpleStruct10;
typedef SimpleStruct<100> SimpleStruct100;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_simple_struct)
{
  BOOST_CHECK(1 == SimpleStruct1().sizeInBytes());
  BOOST_CHECK(sizeInBytes(SimpleStruct1()) == SimpleStruct1().sizeInBytes());
  //
  BOOST_CHECK(10 == SimpleStruct10().sizeInBytes());
  BOOST_CHECK(sizeInBytes(SimpleStruct10()) == SimpleStruct10().sizeInBytes());
  //
  BOOST_CHECK(100 == SimpleStruct100().sizeInBytes());
  BOOST_CHECK(sizeInBytes(SimpleStruct100()) == SimpleStruct100().sizeInBytes());
}

BOOST_AUTO_TEST_CASE(test_std_vector)
{
  std::vector<SimpleStruct1> vector(100);
  BOOST_CHECK(sizeInBytes(vector) == vector.size() * vector[0].sizeInBytes());
}

BOOST_AUTO_TEST_CASE(test_std_array)
{
  std::array<SimpleStruct1, 100> array;
  BOOST_CHECK(sizeInBytes(array) == array.size() * array[0].sizeInBytes());
}

BOOST_AUTO_TEST_SUITE_END()
