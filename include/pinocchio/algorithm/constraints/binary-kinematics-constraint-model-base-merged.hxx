//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_relative_constraint_model_base_hxx__
#define __pinocchio_algorithm_constraints_relative_constraint_model_base_hxx__

namespace pinocchio
{

  template<typename Derived>
  template<int OtherOptions, template<typename, int> class JointCollectionTpl>
  void BinaryKinematicsConstraintModelBase<Derived>::init(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
  {
    nv = model.nv;
    depth_joint1 = static_cast<size_t>(model.supports[joint1_id].size());
    depth_joint2 = static_cast<size_t>(model.supports[joint2_id].size());

    typedef ModelTpl<Scalar, OtherOptions, JointCollectionTpl> Model;
    typedef typename Model::JointModel JointModel;
    static const bool default_sparsity_value = false;
    colwise_joint1_sparsity.fill(default_sparsity_value);
    colwise_joint2_sparsity.fill(default_sparsity_value);

    joint1_span_indexes.reserve(size_t(model.njoints));
    joint2_span_indexes.reserve(size_t(model.njoints));

    JointIndex current1_id = joint1_id > 0 ? joint1_id : JointIndex(0);
    JointIndex current2_id = joint2_id > 0 ? joint2_id : JointIndex(0);

    while (current1_id != current2_id)
    {
      if (current1_id > current2_id)
      {
        const JointModel & joint1 = model.joints[current1_id];
        joint1_span_indexes.push_back((Eigen::Index)current1_id);
        Eigen::Index current1_col_id = joint1.idx_v();
        for (int k = 0; k < joint1.nv(); ++k, ++current1_col_id)
        {
          colwise_joint1_sparsity[current1_col_id] = true;
        }
        current1_id = model.parents[current1_id];
      }
      else
      {
        const JointModel & joint2 = model.joints[current2_id];
        joint2_span_indexes.push_back((Eigen::Index)current2_id);
        Eigen::Index current2_col_id = joint2.idx_v();
        for (int k = 0; k < joint2.nv(); ++k, ++current2_col_id)
        {
          colwise_joint2_sparsity[current2_col_id] = true;
        }
        current2_id = model.parents[current2_id];
      }
    }
    assert(current1_id == current2_id && "current1_id should be equal to current2_id");

    {
      JointIndex current_id = current1_id;
      while (current_id > 0)
      {
        const JointModel & joint = model.joints[current_id];
        joint1_span_indexes.push_back((Eigen::Index)current_id);
        joint2_span_indexes.push_back((Eigen::Index)current_id);
        Eigen::Index current_row_id = joint.idx_v();
        for (int k = 0; k < joint.nv(); ++k, ++current_row_id)
        {
          colwise_joint1_sparsity[current_row_id] = true;
          colwise_joint2_sparsity[current_row_id] = true;
        }
        current_id = model.parents[current_id];
      }
    }
    std::reverse(joint1_span_indexes.begin(), joint1_span_indexes.end());
    std::reverse(joint2_span_indexes.begin(), joint2_span_indexes.end());
    colwise_span_indexes.reserve((size_t)model.nv);
    colwise_sparsity.resize(model.nv);
    colwise_sparsity.setZero();
    for (Eigen::Index col_id = 0; col_id < model.nv; ++col_id)
    {
      if (colwise_joint1_sparsity[col_id] || colwise_joint2_sparsity[col_id])
      {
        colwise_span_indexes.push_back(col_id);
        colwise_sparsity[col_id] = true;
      }
    }

    // Set compliance and Baumgarte parameters.
    m_compliance = ResidualVectorType::Zero(residualSize());
    m_baumgarte_parameters = BaumgarteCorrectorParameters();
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_relative_constraint_model_base_hxx__
