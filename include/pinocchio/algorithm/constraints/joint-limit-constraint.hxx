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
    const JointIndexVector & _activable_joints,
    const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
    const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
    const Eigen::MatrixBase<VectorMarginConfiguration> & margin)
  {
    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;

    PINOCCHIO_CHECK_ARGUMENT_SIZE(lb.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(ub.size(), model.nq);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(margin.size(), model.nq);

    // Check validity of _activable_joints input
    for (const JointIndex joint_id : _activable_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(joint_id > 0, "joint_id should not be equal to zero.");
    }

    //    PINOCCHIO_CHECK_INPUT_ARGUMENT(
    //      check__activable_joints(model, _activable_joints) == -1,
    //      "One of the joint is not supported by JointLimitConstraintModelTpl.")

    // TODO: Should we reserve some activable quantities ?

    // Loop on all q components of activable jointds to identify activable lower and upper
    // constraints, and for each track row_id of related activable joint, idx_q in the configuration
    // and idx_q_reduce in the subpart of q due to activable joints
    VectorOfSize & activable_idx_rows_lower = activable_idx_rows;
    VectorOfSize activable_idx_rows_upper;

    EigenIndexVector & activable_idx_qs_reduce_lower = activable_idx_qs_reduce;
    EigenIndexVector activable_idx_qs_reduce_upper;

    EigenIndexVector & activable_idx_qs_lower = activable_idx_qs;
    EigenIndexVector activable_idx_qs_upper;

    // Prepare the structure to compute sparsity pattern
    EigenIndexVector extended_support;
    extended_support.reserve(size_t(model.nv));

    size_t idx_row =
      0; // there is unique idx_row per joint. In this way, we can have a single row which defines
         // the sparsity pattern associated with the joint, and stored in row_indexes
    nq_reduce = 0;
    for (const JointIndex joint_id : _activable_joints)
    {
      const JointModel & jmodel = model.joints[joint_id];

      const int idx_q = jmodel.idx_q();
      const int idx_v = jmodel.idx_v();
      const int nq = jmodel.nq();
      const int nv = jmodel.nv();
      const auto & has_configuration_limit = jmodel.hasConfigurationLimit();

      bool is_joint_really_active = false; // if at least one of its bound values is finite
      for (int j_qi = 0; j_qi < nq; ++j_qi)
      {
        if (!has_configuration_limit[size_t(j_qi)])
          continue;

        const int q_index = idx_q + j_qi; // index in the plain joint configuration vector q
        const int q_reduce_index = nq_reduce + j_qi;

        if (!(lb[q_index] == -std::numeric_limits<Scalar>::max()
              || lb[q_index] == -std::numeric_limits<Scalar>::infinity()))
        {
          is_joint_really_active = true;
          activable_idx_rows_lower.push_back(idx_row);
          activable_idx_qs_lower.push_back(q_index);
          activable_idx_qs_reduce_lower.push_back(q_reduce_index);
        }
        if (!(ub[q_index] == +std::numeric_limits<Scalar>::max()
              || ub[q_index] == +std::numeric_limits<Scalar>::infinity()))
        {
          is_joint_really_active = true;
          activable_idx_rows_upper.push_back(idx_row);
          activable_idx_qs_upper.push_back(q_index);
          activable_idx_qs_reduce_upper.push_back(q_reduce_index);
        }
      }

      // At least one lower or upper constraint for a component of the joint is active so update the
      // quantity
      if (is_joint_really_active)
      {
        activable_joints.push_back(joint_id);
        idx_row += 1;
        nq_reduce += nq;

        // Compute the sparsity pattern of the joint
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
        row_indexes.push_back(extended_support);
      }
    }

    // Fill row_activable_sparsity_pattern from row_activable_indexes content
    row_sparsity_pattern.resize(row_indexes.size(), BooleanVector::Zero(model.nv));
    for (size_t joint_id = 0; joint_id < row_indexes.size(); ++joint_id)
    {
      auto & sparsity_pattern = row_sparsity_pattern[joint_id];
      const auto & extended_support = row_indexes[joint_id];
      for (const auto val : extended_support)
        sparsity_pattern[val] = true;
    }

    // Recover sizes of constraints
    lower_activable_size = static_cast<int>(activable_idx_rows_lower.size());
    int upper_activable_size = static_cast<int>(activable_idx_rows_upper.size());
    int activable_size = lower_activable_size + upper_activable_size;
    PINOCCHIO_ONLY_USED_FOR_DEBUG(activable_size);

    // Recompose one vectors for all constraint with the convention lower | upper
    activable_idx_rows.insert(
      activable_idx_rows.end(), activable_idx_rows_upper.begin(), activable_idx_rows_upper.end());
    activable_idx_qs_reduce.insert(
      activable_idx_qs_reduce.end(), activable_idx_qs_reduce_upper.begin(),
      activable_idx_qs_reduce_upper.end());
    activable_idx_qs.insert(
      activable_idx_qs.end(), activable_idx_qs_upper.begin(), activable_idx_qs_upper.end());
    assert(size() == activable_size);

    // Fill bound limit and margin for lower and upper constraints
    // Another strategy could be to query the model again but it is not coherent with the existing
    // constructors.
    position_limit = VectorXs::Zero(Eigen::DenseIndex(size()));
    position_margin = VectorXs::Zero(Eigen::DenseIndex(size()));
    Eigen::DenseIndex bound_row_id = 0;
    for (bound_row_id = 0; bound_row_id < lowerSize(); ++bound_row_id)
    {
      const auto activable_idx_q = activable_idx_qs[size_t(bound_row_id)];
      position_limit[bound_row_id] = lb[activable_idx_q];
      assert(margin[activable_idx_q] >= 0);
      position_margin[bound_row_id] = margin[activable_idx_q];
    }
    for (; bound_row_id < size(); ++bound_row_id)
    {
      const auto activable_idx_q = activable_idx_qs[size_t(bound_row_id)];
      position_limit[bound_row_id] = ub[activable_idx_q];
      assert(margin[activable_idx_q] >= 0);
      position_margin[bound_row_id] = margin[activable_idx_q];
    }

    nv_max_atom = 1;
    for (const auto joint_id : activable_joints)
    {
      const auto joint_nv = joint_nvs[joint_id];
      nv_max_atom = std::max(nv_max_atom, joint_nv);
    }

    m_compliance = ComplianceVectorType::Zero(size());
    m_baumgarte_parameters = BaumgarteCorrectorParameters();
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::resize(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & /* model */,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata)
  {
    // Compute notably the constraint constraint_residual
    // This allows to compute which limits are active in the current configuration (data.q_in) which
    // corresponds to the current active set.
    auto & activable_constraint_residual = cdata.activable_constraint_residual;

    auto & active_set_indexes = cdata.active_set_indexes;
    auto & active_idx_rows = cdata.active_idx_rows;
    auto & active_idx_qs_reduce = cdata.active_idx_qs_reduce;
    auto & lower_active_size = cdata.lower_active_size;

    // Reset values
    active_set_indexes.clear();
    active_idx_rows.clear();
    active_idx_qs_reduce.clear();
    lower_active_size = 0;

    // Fill the constraint residual for all activable constraints and detect the active ones.
    // The convention is lower | upper, and negative | positive constraint so:
    // q_l <= q + TMv <= q_up
    // -TMv + (q_l - q) <= 0
    // -TMv + (q_u - q) >= 0

    // We compute all active quanties
    // active_[idx_rows|idx_qs_reduce|nvs|idx_vs] are store
    // but they are not necessary as they are recoverable from active_set_indexes
    // However it implies double referencing in all jacobian methods
    // And for one call of resize/calc, their can be multiple calls to jacobian methods !

    // Lower bounds
    for (std::size_t i = 0; i < static_cast<std::size_t>(lowerSize()); i++)
    {
      const Eigen::DenseIndex i_ = static_cast<Eigen::DenseIndex>(i);
      const Eigen::DenseIndex idx_q = activable_idx_qs[i];
      activable_constraint_residual[i_] = position_limit[i_] - data.q_in[idx_q];
      if (activable_constraint_residual[i_] >= -position_margin[i_])
      {
        active_set_indexes.push_back(i);
        active_idx_rows.push_back(activable_idx_rows[i]);
        active_idx_qs_reduce.push_back(activable_idx_qs_reduce[i]);
        lower_active_size += 1;
      }
    }

    // Upper bounds
    for (std::size_t i = static_cast<std::size_t>(lowerSize());
         i < static_cast<std::size_t>(size()); i++)
    {
      const Eigen::DenseIndex i_ = static_cast<Eigen::DenseIndex>(i);
      const Eigen::DenseIndex idx_q = activable_idx_qs[i];
      activable_constraint_residual[i_] = position_limit[i_] - data.q_in[idx_q];
      if (activable_constraint_residual[i_] <= position_margin[i_])
      {
        active_set_indexes.push_back(i);
        active_idx_rows.push_back(activable_idx_rows[i]);
        active_idx_qs_reduce.push_back(activable_idx_qs_reduce[i]);
      }
    }

    // Resize the constraint residual/compliance storage to the active set size.
    const int active_size = activeSize(cdata);
    cdata.constraint_residual_storage.resize(active_size);

    // Update the active compliance
    cdata.active_compliance_storage.resize(active_size);
    for (int active_row_index = 0; active_row_index < active_size; active_row_index++)
    {
      const auto extended_index =
        Eigen::DenseIndex(active_set_indexes[std::size_t(active_row_index)]);
      cdata.active_compliance[active_row_index] = m_compliance[extended_index];
    }

    // Resize the constraint set so it corresponds to the active set.
    m_set.resize(activeSize(cdata));
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl>
  void JointLimitConstraintModelTpl<Scalar, Options>::calc(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    ConstraintData & cdata) const
  {
    auto & activable_constraint_residual = cdata.activable_constraint_residual;
    auto & constraint_residual = cdata.constraint_residual;

    // TODO: the const cast is only needed because resize touches `m_set`.
    // Introduce set model/data to avoid that.
    const_cast<JointLimitConstraintModelTpl &>(*this).resize(model, data, cdata);
    const std::size_t active_size = static_cast<std::size_t>(this->activeSize(cdata));

    assert(
      constraint_residual.size() == static_cast<int>(active_size)
      && "The active constraint_residual size in constraint data is different from the constraint "
         "model active size.");

    // Fill the constraint residual for all active constraints.
    for (std::size_t active_row_index = 0; active_row_index < active_size; active_row_index++)
    {
      constraint_residual[int(active_row_index)] =
        activable_constraint_residual[int(cdata.active_set_indexes[active_row_index])];
    }

    // Fill the compact tangent map
    auto & compact_tangent_map = cdata.compact_tangent_map;
    pinocchio::compactTangentMap(model, activable_joints, data.q_in, compact_tangent_map);

    auto & rowise_tangent_map = cdata.rowise_tangent_map;
    assert(rowise_tangent_map.size() == size_t(compact_tangent_map.rows()));

    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];
      rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])] =
        compact_tangent_map.row(cdata.active_idx_qs_reduce[constraint_id]).head(constraint_size);
    }
  }

  template<typename Scalar, int Options>
  template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & /*data*/,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
  {
    JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.rows(), this->activeSize(cdata),
      "The input/output Jacobian matrix does not have the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.cols(), model.nv,
      "The input/output Jacobian matrix does not have the right number of cols.");

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;
    jacobian_matrix.setZero();
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];
      const Eigen::Index idx_vs = joint_idx_vs[joint_id];

      jacobian_matrix.row(Eigen::DenseIndex(constraint_id)).segment(idx_vs, constraint_size) =
        -rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])];
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianMatrixProduct(
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
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), activeSize(cdata));
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;
    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id, ++row_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];
      const Eigen::Index idx_vs = joint_idx_vs[joint_id];

      const auto lazy_product_expression =
        -rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])]
        * mat.middleRows(idx_vs, constraint_size);
      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(row_id).noalias() -= lazy_product_expression;
      else
        res.row(row_id).noalias() += lazy_product_expression;
    }
  }

  template<typename Scalar, int Options>
  template<
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointLimitConstraintModelTpl<Scalar, Options>::jacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), activeSize(cdata));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;
    Eigen::DenseIndex row_id = 0;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id, ++row_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];
      const Eigen::Index idx_vs = joint_idx_vs[joint_id];

      const auto lazy_product_expression =
        -rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])].transpose()
        * mat.row(row_id);
      if (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.middleRows(idx_vs, constraint_size).noalias() -= lazy_product_expression;
      else
        res.middleRows(idx_vs, constraint_size).noalias() += lazy_product_expression;
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename VectorNLike,
    ReferenceFrame rf>
  void JointLimitConstraintModelTpl<Scalar, Options>::appendCouplingConstraintInertias(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
    const ReferenceFrameTag<rf> reference_frame) const
  {
    PINOCCHIO_ONLY_USED_FOR_DEBUG(model);
    PINOCCHIO_UNUSED_VARIABLE(reference_frame);

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      diagonal_constraint_inertia.size(), activeSize(cdata),
      "The diagonal_constraint_inertia is of wrong size.");

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];

      const auto & constraint_damping_value =
        diagonal_constraint_inertia[Eigen::DenseIndex(constraint_id)];
      const auto constraint_jacobian =
        -rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])];

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
  void JointLimitConstraintModelTpl<Scalar, Options>::mapConstraintForceToJointTorques(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
    const Eigen::MatrixBase<JointTorquesLike> & joint_torques_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), activeSize(cdata));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_torques_.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);

    auto & joint_torques = joint_torques_.const_cast_derived();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];
      const Eigen::Index idx_vs = joint_idx_vs[joint_id];

      const auto constraint_jacobian =
        -rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])];

      joint_torques.middleRows(idx_vs, constraint_size).noalias() +=
        constraint_jacobian.transpose() * constraint_forces.row(Eigen::DenseIndex(constraint_id));
    }
  }

  template<typename Scalar, int Options>
  template<
    template<typename, int> class JointCollectionTpl,
    typename JointMotionsLike,
    typename ConstraintMotionsLike>
  void JointLimitConstraintModelTpl<Scalar, Options>::mapJointMotionsToConstraintMotion(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
    const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions_.rows(), activeSize(cdata));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);

    auto & constraint_motions = constraint_motions_.const_cast_derived();
    constraint_motions.setZero();

    const auto & rowise_tangent_map = cdata.rowise_tangent_map;
    for (size_t constraint_id = 0; constraint_id < static_cast<std::size_t>(activeSize(cdata));
         ++constraint_id)
    {
      const JointIndex joint_id =
        activable_joints[cdata.active_idx_rows[constraint_id]]; // joint index associated with the
                                                                // constraint
      const Eigen::Index constraint_size = joint_nvs[joint_id];
      const Eigen::Index idx_vs = joint_idx_vs[joint_id];

      const auto constraint_jacobian =
        -rowise_tangent_map[size_t(cdata.active_idx_qs_reduce[constraint_id])];

      constraint_motions.row(Eigen::DenseIndex(constraint_id)).noalias() +=
        constraint_jacobian * joint_motions.middleRows(idx_vs, constraint_size);
    }
  }
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hxx__
