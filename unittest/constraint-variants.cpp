//
// Copyright (c) 2023-2024 INRIA
//

#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include "constraints/init_constraints.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_variants)
{
  Model model;
  buildModels::humanoidRandom(model, true);

  Data data(model);

  PointContactModel rcm = init_constraint<PointContactModel>(model);
  PointContactData rcd(rcm);

  ConstraintModel::ConstraintModelVariant constraint_model_variant = rcm;
  ConstraintModel constraint_model(rcm);
  ConstraintModel constraint_model_equal = rcm;

  ConstraintData constraint_data = rcm.createData();
}

BOOST_AUTO_TEST_CASE(constraint_visitors)
{
  Model model;
  buildModels::humanoidRandom(model, true);

  Data data(model);

  PointContactModel rcm = init_constraint<PointContactModel>(model);
  PointContactData rcd(rcm);

  BOOST_CHECK(ConstraintData(rcd) == ConstraintData(rcd));
  BOOST_CHECK(ConstraintData(rcd) == rcd);

  ConstraintModel constraint_model(rcm);
  ConstraintData constraint_data(rcd);

  // Test size
  {
    BOOST_CHECK(constraint_model.maxResidualSize() == rcm.maxResidualSize());
    BOOST_CHECK(constraint_model.residualSize(constraint_data) == rcm.residualSize(rcd));
  }

  // Test create data visitor
  {
    PointContactData rcd(rcm);
    ConstraintData constraint_data = visitors::createData(constraint_model);
    constraint_data = rcd;
    BOOST_CHECK(constraint_data == rcd);
  }

  // Test calc visitor
  {
    ConstraintData constraint_data1(rcm.createData());
    visitors::calc(constraint_model, model, data, constraint_data1);
    rcm.calc(model, data, rcd);
    //    BOOST_CHECK(rcd == constraint_data1);
    ConstraintData constraint_data2(rcm.createData());
    constraint_model.calc(model, data, constraint_data2);
    //    BOOST_CHECK(rcd == constraint_data2);
  }

  // Test jacobian visitor
  {
    ConstraintData constraint_data(rcm.createData());
    Data::MatrixXs jacobian_matrix1 =
                     Data::MatrixXs::Zero(rcm.residualSize(constraint_data), model.nv),
                   jacobian_matrix2 =
                     Data::MatrixXs::Zero(rcm.residualSize(constraint_data), model.nv),
                   jacobian_matrix_ref =
                     Data::MatrixXs::Zero(rcm.residualSize(constraint_data), model.nv);
    rcm.jacobian(model, data, rcd, jacobian_matrix_ref);
    visitors::jacobian(constraint_model, model, data, constraint_data, jacobian_matrix1);
    BOOST_CHECK(jacobian_matrix1 == jacobian_matrix_ref);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix2);
    BOOST_CHECK(jacobian_matrix2 == jacobian_matrix_ref);
  }

  // Test getRowIndexes
  {
    ConstraintData constraint_data(rcm.createData());
    rcm.calc(model, data, rcd);
    for (Eigen::Index row_id = 0; row_id < constraint_model.residualSize(rcd); ++row_id)
    {
      BOOST_CHECK(
        constraint_model.getRowIndexes(constraint_data, row_id) == rcm.getRowIndexes(rcd, row_id));
    }
  }

  // Test getRowSparsityPattern
  {
    ConstraintData constraint_data(rcm.createData());
    for (Eigen::Index row_id = 0; row_id < constraint_model.residualSize(rcd); ++row_id)
    {
      BOOST_CHECK(
        constraint_model.getRowSparsityPattern(constraint_data, row_id)
        == rcm.getRowSparsityPattern(rcd, row_id));
    }
  }

  // Test jacobianMatrixProduct
  {
    const Eigen::Index num_cols = 20;
    ConstraintData constraint_data(rcm.createData());
    const Data::MatrixXs input_matrix = Data::MatrixXs::Random(model.nv, num_cols);
    Data::MatrixXs output_matrix1(rcm.residualSize(constraint_data), num_cols),
      output_matrix2(rcm.residualSize(constraint_data), num_cols),
      output_matrix_ref(rcm.residualSize(constraint_data), num_cols);
    rcm.jacobianMatrixProduct(model, data, rcd, input_matrix, output_matrix_ref);
    visitors::jacobianMatrixProduct(
      constraint_model, model, data, constraint_data, input_matrix, output_matrix1);
    BOOST_CHECK(output_matrix1 == output_matrix_ref);
    constraint_model.jacobianMatrixProduct(
      model, data, constraint_data, input_matrix, output_matrix2);
    BOOST_CHECK(output_matrix2 == output_matrix_ref);
  }

  // Test jacobianTransposeMatrixProduct
  {
    const Eigen::Index num_cols = 20;
    ConstraintData constraint_data(rcm.createData());
    const Data::MatrixXs input_matrix =
      Data::MatrixXs::Random(rcm.residualSize(constraint_data), num_cols);
    Data::MatrixXs output_matrix1(model.nv, num_cols), output_matrix2(model.nv, num_cols),
      output_matrix_ref(model.nv, num_cols);
    rcm.jacobianTransposeMatrixProduct(model, data, rcd, input_matrix, output_matrix_ref);
    visitors::jacobianTransposeMatrixProduct(
      constraint_model, model, data, constraint_data, input_matrix, output_matrix1);
    BOOST_CHECK(output_matrix1 == output_matrix_ref);
    constraint_model.jacobianTransposeMatrixProduct(
      model, data, constraint_data, input_matrix, output_matrix2);
    BOOST_CHECK(output_matrix2 == output_matrix_ref);
  }
}

BOOST_AUTO_TEST_SUITE_END()
