//
// Copyright (c) 2024-2025 INRIA
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/constraints/joint-friction-constraint.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include <iostream>

// Helpers
#include "constraints/jacobians-checker.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;
typedef JointFrictionConstraintModel::EigenIndexVector EigenIndexVector;
typedef JointFrictionConstraintModel::BooleanVector BooleanVector;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_empty_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Data data(model);

  const Model::IndexVector empty_active_joint_ids;

  JointFrictionConstraintModel constraint(model, empty_active_joint_ids);
}

BOOST_AUTO_TEST_CASE(constraint_constructor)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Data data(model);
  const auto & parents_fromRow = data.parents_fromRow;

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  JointFrictionConstraintModel constraint(model, active_joint_ids);
  JointFrictionConstraintData constraint_data = constraint.createData();

  // Check size
  {
    int total_size = 0;
    for (const JointIndex joint_id : active_joint_ids)
    {
      total_size += model.joints[joint_id].nv();
    }
    BOOST_CHECK(constraint.maxResidualSize() == total_size);
    BOOST_CHECK(constraint.getActiveDofs().size() == size_t(total_size));
  }

  // Check sparsity pattern
  {
    const EigenIndexVector & active_dofs = constraint.getActiveDofs();
    for (size_t row_id = 0; row_id < size_t(constraint.maxResidualSize()); ++row_id)
    {
      const Eigen::Index dof_id = active_dofs[row_id];
      const BooleanVector & row_sparsity_pattern =
        constraint.getRowSparsityPattern(model, data, constraint_data, Eigen::Index(row_id));
      const EigenIndexVector & row_active_indexes =
        constraint.getRowIndexes(model, data, constraint_data, Eigen::Index(row_id));

      // Check that the rest of the indexes greater than dof_id are not active.
      BOOST_CHECK((row_sparsity_pattern.tail(model.nv - 1 - dof_id).array() == false).all());

      Eigen::Index id = dof_id;
      while (parents_fromRow[size_t(id)] > -1)
      {
        BOOST_CHECK(row_sparsity_pattern[id] == true);
        id = parents_fromRow[size_t(id)];
      }

      for (const Eigen::Index active_id : row_active_indexes)
      {
        BOOST_CHECK(row_sparsity_pattern[active_id] == true);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(cast)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Eigen::VectorXd q = neutral(model);

  const Data data(model);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  JointFrictionConstraintModel cm(model, active_joint_ids);

  const auto cm_cast_double = cm.cast<double>();
  BOOST_CHECK(cm_cast_double == cm);

  const auto cm_cast_long_double = cm.cast<long double>();
  BOOST_CHECK(cm_cast_long_double.cast<double>() == cm);
}

BOOST_AUTO_TEST_CASE(constraint_jacobian)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Eigen::VectorXd q = neutral(model);

  const Data data(model);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  JointFrictionConstraintModel constraint_model(model, active_joint_ids);
  JointFrictionConstraintData constraint_data(constraint_model);

  Eigen::MatrixXd jacobian_matrix(constraint_model.maxResidualSize(), model.nv);
  constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);

  const EigenIndexVector & active_dofs = constraint_model.getActiveDofs();
  for (Eigen::Index row_id = 0; row_id < constraint_model.maxResidualSize(); ++row_id)
  {
    const Eigen::Index dof_id = active_dofs[size_t(row_id)];
    BOOST_CHECK(jacobian_matrix.row(row_id).sum() == 1.);
    BOOST_CHECK(jacobian_matrix(row_id, dof_id) == 1.);
    BOOST_CHECK(
      (dof_id - 1) > 0 ? (jacobian_matrix.row(row_id).head(dof_id - 1).array() == 0).all() : true);
    BOOST_CHECK(
      (model.nv - dof_id - 1) > 0
        ? (jacobian_matrix.row(row_id).tail(model.nv - dof_id - 1).array() == 0).all()
        : true);
  }

  check_jacobians_operations(model, data, constraint_model, constraint_data);
}

BOOST_AUTO_TEST_CASE(constraint_coupling_inertia)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);

  const Eigen::VectorXd q = neutral(model);

  Data data(model);
  computeJointJacobians(model, data, q);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  JointFrictionConstraintModel constraint_model(model, active_joint_ids);
  JointFrictionConstraintData constraint_data(constraint_model);

  constraint_model.calc(model, data, constraint_data);
  const Eigen::VectorXd diagonal_inertia =
    Eigen::VectorXd::Random(constraint_model.residualSize(constraint_data)).array().square();
  constraint_model.appendCouplingConstraintInertias(
    model, data, constraint_data, diagonal_inertia, WorldFrameTag());

  Eigen::Index row_id = 0;

  for (const auto joint_id : active_joint_ids)
  {
    //    std::cout << "joint_id: " << joint_id << std::endl;

    const auto & jmodel = model.joints[joint_id];
    const auto jmodel_nv = jmodel.nv();
    // const auto jmodel_idx_v = jmodel.idx_v();

    const auto diagonal_inertia_segment = diagonal_inertia.segment(row_id, jmodel_nv);

    BOOST_CHECK(diagonal_inertia_segment == data.joint_apparent_inertia[joint_id].diagonal());

    row_id += jmodel_nv;
    //    std::cout << "----" << std::endl;
  }

  Eigen::MatrixXd jacobian_matrix(constraint_model.residualSize(constraint_data), model.nv);
  constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);

  const Eigen::MatrixXd joint_space_constraint_inertia =
    jacobian_matrix.transpose() * diagonal_inertia.asDiagonal() * jacobian_matrix;

  for (const auto joint_id : active_joint_ids)
  {
    const auto & jmodel = model.joints[joint_id];
    const auto jmodel_nv = jmodel.nv();
    const auto jmodel_idx_v = jmodel.idx_v();

    BOOST_CHECK(
      joint_space_constraint_inertia.block(jmodel_idx_v, jmodel_idx_v, jmodel_nv, jmodel_nv)
        .isApprox(data.joint_apparent_inertia[joint_id]));
  }
}

BOOST_AUTO_TEST_CASE(check_maps)
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model, true);
  Data data(model), data_ref(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);

  const std::string RF_name = "rleg6_joint";
  const JointIndex RF_id = model.getJointId(RF_name);

  const Model::IndexVector & RF_support = model.supports[RF_id];
  const Model::IndexVector active_joint_ids(RF_support.begin() + 1, RF_support.end());

  JointFrictionConstraintModel constraint_model(model, active_joint_ids);
  JointFrictionConstraintData constraint_data(constraint_model),
    constraint_data_ref(constraint_model);

  const Eigen::VectorXd q = neutral(model);
  computeJointJacobians(model, data, q);
  constraint_model.calc(model, data, constraint_data);
  computeJointJacobians(model, data_ref, q);
  constraint_model.calc(model, data_ref, constraint_data_ref);

  const auto constraint_jacobian_ref =
    constraint_model.jacobian(model, data_ref, constraint_data_ref);

  // Test mapConstraintForcesToJointTorques
  {
    const Eigen::VectorXd constraint_forces =
      Eigen::VectorXd::Random(constraint_model.residualSize(constraint_data));

    Eigen::VectorXd joint_torques_ref = Eigen::VectorXd::Zero(model.nv);
    joint_torques_ref = constraint_jacobian_ref.transpose() * constraint_forces;

    Eigen::VectorXd joint_torques_ref2 = Eigen::VectorXd::Zero(model.nv);
    constraint_model.jacobianTransposeMatrixProduct(
      model, data_ref, constraint_data_ref, constraint_forces, joint_torques_ref2, SetTo());

    Eigen::VectorXd joint_torques = Eigen::VectorXd::Zero(model.nv);
    constraint_model.mapConstraintForceToJointTorques(
      model, data_ref, constraint_data, constraint_forces, joint_torques);

    BOOST_CHECK(joint_torques.isApprox(joint_torques_ref));
    BOOST_CHECK(joint_torques.isApprox(joint_torques_ref2));
  }

  // Test mapJointMotionsToConstraintMotions
  {
    const Eigen::VectorXd joint_motions = Eigen::VectorXd::Random(model.nv);

    Eigen::VectorXd constraint_motions_ref =
      Eigen::VectorXd::Zero(constraint_model.residualSize(constraint_data));
    constraint_motions_ref = constraint_jacobian_ref * joint_motions;

    Eigen::VectorXd constraint_motions_ref2 =
      Eigen::VectorXd::Zero(constraint_model.residualSize(constraint_data));
    constraint_model.jacobianMatrixProduct(
      model, data_ref, constraint_data_ref, joint_motions, constraint_motions_ref2, SetTo());

    Eigen::VectorXd constraint_motions =
      -Eigen::VectorXd::Ones(constraint_model.residualSize(constraint_data));
    constraint_model.mapJointMotionsToConstraintMotion(
      model, data_ref, constraint_data, joint_motions, constraint_motions);

    BOOST_CHECK(constraint_motions.isApprox(constraint_motions_ref));
    BOOST_CHECK(constraint_motions.isApprox(constraint_motions_ref2));
  }
}

BOOST_AUTO_TEST_SUITE_END()
