//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__

namespace pinocchio
{

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorLowerConfiguration,
    typename VectorUpperConfiguration,
    typename VectorMarginConfiguration>
  void JointLimitConstraintModelTpl<Scalar, Options>::init(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const JointIndexVector & activable_joints,
    const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
    const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
    const Eigen::MatrixBase<VectorMarginConfiguration> & margin)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;

    PINOCCHIO_CHECK_ARGUMENT_SIZE(lb.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(ub.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(margin.size(), model.nq);

    // Check validity of activable_joints input
    for (const JointIndex joint_id : activable_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(joint_id > 0, "joint_id should not be equal to zero.");
    }

    // Loop on all q components of activable jointds to identify activable lower and upper
    // constraints, and for each track row_id of related activable joint, idx_q in the configuration
    // and idx_q_reduce in the subpart of q due to activable joints
    VectorOfSize & activable_idx_in_selected_lower = activable_idx_in_selected;
    VectorOfSize activable_idx_in_selected_upper;

    EigenIndexVector & activable_idx_qs_reduce_lower = activable_idx_qs_reduce;
    EigenIndexVector activable_idx_qs_reduce_upper;

    EigenIndexVector & activable_idx_qs_lower = activable_idx_qs;
    EigenIndexVector activable_idx_qs_upper;

    // Prepare the structure to compute sparsity pattern
    EigenIndexVector extended_support;
    extended_support.reserve(size_t(model.nv));

    size_t idx_selected =
      0; // there is unique idx_row per joint. In this way, we can have a single row which defines
         // the sparsity pattern associated with the joint, and stored in row_indexes
    // Note: Instead of looping on joint, we can use a visitor to loop on atomic joints
    // It correspond to sub-looping on joint composite
    // It was done in the past, but for coherence with the algorithms, that consider the
    // joint composite as a block it has been removed.
    nq_reduce = 0;
    max_of_nvs = 1;
    for (const JointIndex joint_id : activable_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];

      const int idx_q = jmodel.idx_q();
      const int idx_v = jmodel.idx_v();
      const int nq = jmodel.nq();
      const int nv = jmodel.nv();
      const auto & has_configuration_limit = jmodel.hasConfigurationLimit();

      bool is_joint_selected = false; // if at least one of its bound values is finite
      for (int j_qi = 0; j_qi < nq; ++j_qi)
      {
        if (!has_configuration_limit[static_cast<size_t>(j_qi)])
          continue;

        const int q_index = idx_q + j_qi; // index in the plain joint configuration vector q
        const int q_reduce_index = nq_reduce + j_qi;

        if (!(lb[q_index] == -std::numeric_limits<Scalar>::max()
              || lb[q_index] == -std::numeric_limits<Scalar>::infinity()))
        {
          is_joint_selected = true;
          activable_idx_in_selected_lower.push_back(idx_selected);
          activable_idx_qs_reduce_lower.push_back(q_reduce_index);
          activable_idx_qs_lower.push_back(q_index);
        }
        if (!(ub[q_index] == +std::numeric_limits<Scalar>::max()
              || ub[q_index] == +std::numeric_limits<Scalar>::infinity()))
        {
          is_joint_selected = true;
          activable_idx_in_selected_upper.push_back(idx_selected);
          activable_idx_qs_reduce_upper.push_back(q_reduce_index);
          activable_idx_qs_upper.push_back(q_index);
        }
      }

      // At least one lower or upper constraint for a component of the joint is selected
      // so calculate its sparsity pattern
      if (is_joint_selected)
      {
        selected_joints.push_back(joint_id);
        selected_joint_nqs.push_back(nq);
        selected_joint_nvs.push_back(nv);
        selected_joint_idx_vs.push_back(idx_v);

        idx_selected += 1;
        nq_reduce += nq;
        max_of_nvs = std::max(max_of_nvs, nv);

        // Compute the row indexes of the joint
        const auto & joint_support = model.supports[joint_id];
        extended_support.clear();
        for (size_t i = 1; i < joint_support.size() - 1; ++i)
        {
          const JointIndex joint_support_id = joint_support[i];
          const JointModel & joint_support = model.joints[joint_support_id];
          const int joint_support_nv = joint_support.nv();
          const int joint_support_idx_v = joint_support.idx_v();
          for (int k = 0; k < joint_support_nv; ++k)
          {
            const int extended_row_id = joint_support_idx_v + k;
            extended_support.push_back(extended_row_id);
          }
        }
        for (int k = 0; k < nv; ++k)
        {
          const int extended_row_id = idx_v + k;
          extended_support.push_back(extended_row_id);
        }
        selected_row_indexes.push_back(extended_support);
      }
    }

    // Fill selected_row_sparsity_pattern from selected_row_indexes content
    selected_row_sparsity_pattern.resize(
      selected_row_indexes.size(), BooleanVector::Zero(model.nv));
    for (size_t idx_sel = 0; idx_sel < selected_row_indexes.size(); ++idx_sel)
    {
      auto & sparsity_pattern = selected_row_sparsity_pattern[idx_sel];
      const auto & extended_support = selected_row_indexes[idx_sel];
      for (const auto val : extended_support)
        sparsity_pattern[val] = true;
    }

    // Recover max sizes of constraint
    lower_max_residual_size = static_cast<int>(activable_idx_in_selected_lower.size());

    const int upper_max_residual_size = static_cast<int>(activable_idx_in_selected_upper.size());
    PINOCCHIO_ONLY_USED_FOR_DEBUG(upper_max_residual_size);
    const int max_residual_size = lower_max_residual_size + upper_max_residual_size;
    PINOCCHIO_ONLY_USED_FOR_DEBUG(max_residual_size);

    // Recompose one vectors for all constraint with the convention lower | upper
    activable_idx_in_selected.insert(
      activable_idx_in_selected.end(), activable_idx_in_selected_upper.begin(),
      activable_idx_in_selected_upper.end());
    activable_idx_qs_reduce.insert(
      activable_idx_qs_reduce.end(), activable_idx_qs_reduce_upper.begin(),
      activable_idx_qs_reduce_upper.end());
    activable_idx_qs.insert(
      activable_idx_qs.end(), activable_idx_qs_upper.begin(), activable_idx_qs_upper.end());

    assert(maxResidualSize() == max_residual_size);

    // Set activable_[position_limit|margin] of size maxResidualSize from lb, ub, margin of size
    // model.nq
    setPositionLimitAndMargin(lb, ub, margin);

    // Data member
    m_compliance = ComplianceVectorType::Zero(maxResidualSize());
    m_baumgarte_parameters = BaumgarteCorrectorParameters();
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::computeResidualAndActiveConstraints(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & /* model */,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata) const
  {
    // Compute notably the constraint constraint_residual
    // This allows to compute which limits are active in the current configuration (data.q_in) which
    // corresponds to the constraints to consider.

    auto & active_idx_in_activable = cdata.active_idx_in_activable;
    auto & active_idx_in_selected = cdata.active_idx_in_selected;
    auto & active_idx_qs_reduce = cdata.active_idx_qs_reduce;
    auto & lower_residual_size = cdata.lower_residual_size;

    // Reset values
    active_idx_in_activable.clear();
    active_idx_in_selected.clear();
    active_idx_qs_reduce.clear();
    lower_residual_size = 0;

    // Fill the constraint residual for all activable constraints and detect the active ones.
    // The convention is lower | upper, and negative | positive constraint so:
    // q_l <= q + TMv <= q_up
    // TMv + (q - q_l) >= 0 --> J = TM, res = q-q_l
    // -TMv + (q_u - q) >= 0 --> J = -TM, res = q_u-q

    // We track the active constraints among the activable using active_idx_in_activable.
    // We also compute the proxy active_[idx_in_selected|idx_qs_reduce] which are:
    // active_X[i] = activable_X[active_idx_in_activable[i]]
    // to avoid double referencing.

    auto & activable_constraint_residual = cdata.activable_constraint_residual;
    std::size_t i = 0;
    // Lower bounds
    for (; i < static_cast<std::size_t>(lowerMaxResidualSize()); i++)
    {
      const Eigen::Index i_ = static_cast<Eigen::Index>(i);
      const Eigen::Index idx_q = activable_idx_qs[i];
      activable_constraint_residual[i_] = data.q_in[idx_q] - activable_position_limit[i_];
      if (activable_constraint_residual[i_] <= activable_position_margin[i_])
      {
        active_idx_in_activable.push_back(i);
        // Update proxis as well
        active_idx_in_selected.push_back(activable_idx_in_selected[i]);
        active_idx_qs_reduce.push_back(activable_idx_qs_reduce[i]);
        lower_residual_size += 1;
      }
    }
    // Upper bounds
    for (; i < static_cast<std::size_t>(maxResidualSize()); i++)
    {
      const Eigen::Index i_ = static_cast<Eigen::Index>(i);
      const Eigen::Index idx_q = activable_idx_qs[i];
      activable_constraint_residual[i_] = activable_position_limit[i_] - data.q_in[idx_q];
      if (activable_constraint_residual[i_] <= activable_position_margin[i_])
      {
        active_idx_in_activable.push_back(i);
        // Update proxis as well
        active_idx_in_selected.push_back(activable_idx_in_selected[i]);
        active_idx_qs_reduce.push_back(activable_idx_qs_reduce[i]);
      }
    }

    // Resize the constraint residual/compliance storage to the constraint size.
    cdata.constraint_residual_storage.resize(residualSize(cdata));
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::calcImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata) const
  {
    // Compute the residual for all activable constraints and identify the active constraints
    computeResidualAndActiveConstraints(model, data, cdata);
    auto & activable_constraint_residual = cdata.activable_constraint_residual;

    // Fill the compact tangent map for the system in configuration q_in
    auto & compact_tangent_map = cdata.compact_tangent_map;
    pinocchio::compactTangentMap(model, selected_joints, data.q_in, compact_tangent_map);

    // For each constraint, recover the residual and store the rowise tangent map
    const std::size_t csize = static_cast<std::size_t>(residualSize(cdata));
    auto & constraint_residual = cdata.constraint_residual;
    auto & rowise_tangent_map = cdata.rowise_tangent_map;
    assert(
      constraint_residual.size() == int(csize)
      && "The active constraint_residual size in constraint data is different from the constraint "
         "model active size.");
    for (std::size_t constraint_id = 0; constraint_id < csize; constraint_id++)
    {
      const Eigen::Index idx_q_reduce = cdata.active_idx_qs_reduce[constraint_id];
      const Eigen::Index idx_in_activable =
        Eigen::Index(cdata.active_idx_in_activable[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];

      // Recover the constraint residual
      constraint_residual[Eigen::Index(constraint_id)] =
        activable_constraint_residual[idx_in_activable];

      // Store the rowise_tangent_map
      rowise_tangent_map[static_cast<std::size_t>(idx_q_reduce)] =
        compact_tangent_map.row(idx_q_reduce).head(constraint_size);
    }
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & /*data*/,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
  {
    JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.rows(), this->residualSize(cdata),
      "The input/output Jacobian matrix does not have the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.cols(), model.nv,
      "The input/output Jacobian matrix does not have the right number of cols.");

    jacobian_matrix.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;

    std::size_t constraint_id = 0;
    // Lower bounds
    for (; constraint_id < static_cast<size_t>(lowerResidualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      jacobian_matrix.row(Eigen::Index(constraint_id)).segment(idx_v, constraint_size) =
        rowise_tangent_map[idx_q_reduce];
    }
    // Upper bounds
    for (; constraint_id < static_cast<size_t>(residualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      jacobian_matrix.row(Eigen::Index(constraint_id)).segment(idx_v, constraint_size) =
        -rowise_tangent_map[idx_q_reduce];
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianMatrixProductImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), residualSize(cdata));
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;

    std::size_t constraint_id = 0;
    // Lower bounds
    for (; constraint_id < static_cast<size_t>(lowerResidualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto lazy_product_expression =
        rowise_tangent_map[idx_q_reduce] * mat.middleRows(idx_v, constraint_size);
      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(Eigen::Index(constraint_id)).noalias() -= lazy_product_expression;
      else
        res.row(Eigen::Index(constraint_id)).noalias() += lazy_product_expression;
    }
    // Upper bounds
    for (; constraint_id < static_cast<size_t>(residualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto lazy_product_expression =
        -rowise_tangent_map[idx_q_reduce] * mat.middleRows(idx_v, constraint_size);
      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(Eigen::Index(constraint_id)).noalias() -= lazy_product_expression;
      else
        res.row(Eigen::Index(constraint_id)).noalias() += lazy_product_expression;
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianTransposeMatrixProductImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), residualSize(cdata));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;

    std::size_t constraint_id = 0;
    // Lower bounds
    for (; constraint_id < static_cast<size_t>(lowerResidualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto lazy_product_expression =
        rowise_tangent_map[idx_q_reduce].transpose() * mat.row(Eigen::Index(constraint_id));
      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.middleRows(idx_v, constraint_size).noalias() -= lazy_product_expression;
      else
        res.middleRows(idx_v, constraint_size).noalias() += lazy_product_expression;
    }
    // Upper bounds
    for (; constraint_id < static_cast<size_t>(residualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto lazy_product_expression =
        -rowise_tangent_map[idx_q_reduce].transpose() * mat.row(Eigen::Index(constraint_id));
      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.middleRows(idx_v, constraint_size).noalias() -= lazy_product_expression;
      else
        res.middleRows(idx_v, constraint_size).noalias() += lazy_product_expression;
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorNLike,
    ReferenceFrame rf>
  void JointLimitConstraintModelTpl<Scalar, Options>::appendCouplingConstraintInertiasImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
    const ReferenceFrameTag<rf> reference_frame) const
  {
    PINOCCHIO_ONLY_USED_FOR_DEBUG(model);
    PINOCCHIO_UNUSED_VARIABLE(reference_frame);

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      diagonal_constraint_inertia.size(), residualSize(cdata),
      "The diagonal_constraint_inertia is of wrong size.");

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;

    // Lower bounds and upper bounds together as (-R)^T(-R) = R^TR
    for (std::size_t constraint_id = 0; constraint_id < static_cast<size_t>(residualSize(cdata));
         constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const JointIndex joint_id = selected_joints[idx_in_selected];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];

      const auto & constraint_damping_value =
        diagonal_constraint_inertia[Eigen::Index(constraint_id)];
      const auto constraint_jacobian = rowise_tangent_map[idx_q_reduce];

      assert(
        joint_id > 0 && joint_id < JointIndex(model.njoints) && "joint_id value is incorrect.");

      auto & support_joint_apparent_inertia = data.joint_apparent_inertia[joint_id];
      PINOCCHIO_ONLY_USED_FOR_DEBUG(constraint_size);
      assert(support_joint_apparent_inertia.rows() == constraint_size);
      assert(support_joint_apparent_inertia.cols() == constraint_size);

      support_joint_apparent_inertia.noalias() +=
        constraint_damping_value * constraint_jacobian.transpose() * constraint_jacobian;
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename ConstraintForcesLike,
    typename JointTorquesLike>
  void JointLimitConstraintModelTpl<Scalar, Options>::mapConstraintForceToJointTorquesImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
    const Eigen::MatrixBase<JointTorquesLike> & joint_torques_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), residualSize(cdata));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_torques_.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);

    auto & joint_torques = joint_torques_.const_cast_derived();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;

    std::size_t constraint_id = 0;
    // Lower bounds
    for (; constraint_id < static_cast<std::size_t>(lowerResidualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto constraint_jacobian = rowise_tangent_map[idx_q_reduce];

      joint_torques.middleRows(idx_v, constraint_size).noalias() +=
        constraint_jacobian.transpose() * constraint_forces.row(Eigen::Index(constraint_id));
    }
    // Upper bounds
    for (; constraint_id < static_cast<std::size_t>(residualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto constraint_jacobian = -rowise_tangent_map[idx_q_reduce];

      joint_torques.middleRows(idx_v, constraint_size).noalias() +=
        constraint_jacobian.transpose() * constraint_forces.row(Eigen::Index(constraint_id));
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename JointMotionsLike,
    typename ConstraintMotionsLike>
  void JointLimitConstraintModelTpl<Scalar, Options>::mapJointMotionsToConstraintMotionImpl(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
    const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions_.rows(), residualSize(cdata));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);

    auto & constraint_motions = constraint_motions_.const_cast_derived();
    constraint_motions.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;

    std::size_t constraint_id = 0;
    // Lower bounds
    for (; constraint_id < static_cast<std::size_t>(lowerResidualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto constraint_jacobian = rowise_tangent_map[idx_q_reduce];

      constraint_motions.row(Eigen::Index(constraint_id)).noalias() +=
        constraint_jacobian * joint_motions.middleRows(idx_v, constraint_size);
    }
    // Upper bounds
    for (; constraint_id < static_cast<std::size_t>(residualSize(cdata)); constraint_id++)
    {
      const size_t idx_q_reduce = static_cast<size_t>(cdata.active_idx_qs_reduce[constraint_id]);
      const size_t idx_in_selected = cdata.active_idx_in_selected[constraint_id];
      const Eigen::Index constraint_size = selected_joint_nvs[idx_in_selected];
      const Eigen::Index idx_v = selected_joint_idx_vs[idx_in_selected];

      const auto constraint_jacobian = -rowise_tangent_map[idx_q_reduce];

      constraint_motions.row(Eigen::Index(constraint_id)).noalias() +=
        constraint_jacobian * joint_motions.middleRows(idx_v, constraint_size);
    }
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
