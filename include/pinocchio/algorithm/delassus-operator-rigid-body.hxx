//
// Copyright (c) 2024-2026 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/algorithm/constraints/constraints.hpp"

#include "pinocchio/algorithm/loop-constrained-aba.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/constraints/utils.hpp"

#include "pinocchio/algorithm/delassus-operator-rigid-body-visitors.hxx"

namespace pinocchio
{

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::
    rebuild(
      const ModelHolder & model_ref,
      const DataHolder & data_ref,
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref)
  {
    // Rebuild quantities related to model
    m_model_ref = model_ref;
    m_data_ref = data_ref;
    m_internal_data.rebuild(model());

    // Rebuild quantities related to constraints
    m_constraint_models_ref = constraint_models_ref;
    m_constraint_datas_ref = constraint_datas_ref;

    m_size = residualSize(helper::get_ref(constraint_models_ref));

    // resize quantities
    m_damping = VectorXs::Constant(m_size, m_min_damping_value).asDiagonal();
    m_compliance_storage.resize(m_size);
    m_compliance.setZero();
    m_sum_compliance_damping = VectorXs::Constant(m_size, m_min_damping_value).asDiagonal();
    m_sum_compliance_damping_inverse = VectorXs::Zero(m_size).asDiagonal();

    assert(m_damping.rows() == m_size);
    assert(m_damping.cols() == m_size);
    assert(m_compliance.size() == m_size);
    assert(m_sum_compliance_damping.rows() == m_size);
    assert(m_sum_compliance_damping_inverse.rows() == m_size);

    retrieveConstraintCompliance(helper::get_ref(constraint_models_ref), m_compliance);

    computeJointMinimalOrdering(model(), data(), helper::get_ref(constraint_models_ref));
    updateSumComplianceDamping();
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::compute_or_update_decomposition(bool apply_on_the_right, bool solve_in_place)
  {
    typedef typename Data::Inertia Inertia;
    using Matrix6 = typename Inertia::Matrix6;

    const Model & model_ref = model();
    Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();

    for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
    {
      const auto & joint_inertia = model_ref.inertias[i];
      if (apply_on_the_right)
        data_ref.Yaba[i] = joint_inertia.matrix();
      if (solve_in_place)
      {
        const Inertia oinertia = data_ref.oMi[i].act(joint_inertia);
        data_ref.oYaba_augmented[i] = oinertia.matrix();
      }
    }

    if (solve_in_place)
    {
      for (JointIndex joint_id = 1; joint_id < JointIndex(model_ref.njoints); ++joint_id)
      {
        const auto joint_nv = model_ref.nvs[joint_id];
        const auto joint_idx_v = model_ref.idx_vs[joint_id];

        data_ref.joint_apparent_inertia[joint_id] =
          model_ref.armature.segment(joint_idx_v, joint_nv).asDiagonal();
      }

      data_ref.joint_cross_coupling.apply([](Matrix6 & v) { v.setZero(); });

      // Append constraint inertia to oYaba_augmented
      m_sum_compliance_damping_inverse = m_sum_compliance_damping.inverse();
      assert(!m_sum_compliance_damping_inverse.hasNaN());

      const auto & blocks = m_sum_compliance_damping_inverse.blocks();
      if (blocks.size() == 1 && constraint_models_ref.size() > 1) // we assume we have a single
                                                                  // diagonal block to dispatch on
                                                                  // all the contraints
      {
        const auto & diagonal_block = blocks[0];
        assert(diagonal_block.type() == MatrixBlockType::Diagonal);
        Eigen::Index row_id = 0;
        typedef typename BlockDiagonalMatrix::ConstVectorMap ConstVectorMap;
        const auto & compliance_damping_inverse_vector =
          remap<ConstVectorMap>(diagonal_block.container());
        for (std::size_t constraint_id = 0; constraint_id < constraint_models_ref.size();
             ++constraint_id)
        {
          const auto & cmodel = helper::get_ref(constraint_models_ref[constraint_id]);
          const auto & cdata = helper::get_ref(constraint_datas_ref[constraint_id]);

          const auto constraint_size = cmodel.residualSize();
          const auto constraint_diagonal_inertia =
            compliance_damping_inverse_vector.segment(row_id, constraint_size);

          cmodel.appendCouplingConstraintInertias(
            model_ref, data_ref, cdata, constraint_diagonal_inertia, WorldFrameTag());
          row_id += constraint_size;
        }
        assert(row_id == size());
      }
      else // we have block diagonal matrix, each block being assigned to a constraint
      {
        std::size_t inner_constraint_id = 0;
        for (std::size_t constraint_id = 0; constraint_id < constraint_models_ref.size();
             ++constraint_id)
        {
          const auto & cmodel = helper::get_ref(constraint_models_ref[constraint_id]);
          const auto & cdata = helper::get_ref(constraint_datas_ref[constraint_id]);
          cmodel.appendCouplingConstraintInertias(
            model_ref, data_ref, cdata, blocks, WorldFrameTag(), inner_constraint_id);
        }
      }
    }

#define DO_PASS(apply_on_the_right_v, solve_in_place_v)                                            \
  {                                                                                                \
    typedef DelassusOperatorRigidBodySystemsComputeBackwardPass<                                   \
      DelassusOperatorRigidBodySystemsTpl, apply_on_the_right_v, solve_in_place_v>                 \
      Pass2;                                                                                       \
    for (const JointIndex i : data_ref.joint_elimination_order)                                    \
    {                                                                                              \
      typename Pass2::ArgsType args(model_ref, data_ref);                                          \
      Pass2::run(model_ref.joints[i], data_ref.joints[i], args);                                   \
    }                                                                                              \
  }

    if (apply_on_the_right)
    {
      if (solve_in_place)
      {
        DO_PASS(true, true);
      }
      else
      {
        DO_PASS(true, false);
      }
    }
    else
    {
      if (solve_in_place)
      {
        DO_PASS(false, true);
      }
      else
      {
        DO_PASS(false, false);
      }
    }
#undef DO_PASS

    if (solve_in_place)
      m_solve_in_place_dirty = false;
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  template<typename MatrixIn, typename MatrixOut>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::
    applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & rhs,
      const Eigen::MatrixBase<MatrixOut> & res_,
      bool with_damping) const
  {
    MatrixOut & res = res_.const_cast_derived();
    PINOCCHIO_CHECK_SAME_MATRIX_SIZE(rhs, res);

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();
    auto & internal_data = this->m_internal_data;
    auto & u = internal_data.u;

    // Make a pass over the whole set of constraints to add the contributions of constraint forces
    // u and internal_data.f are reset by mapConstraintForcesToJointSpace
    mapConstraintForcesToJointSpace(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, rhs, m_internal_data.f, u,
      LocalFrameTag());
    // TODO(jcarpent): extend the code to operator on matrices

    //    typedef Eigen::Map<VectorXs,EIGEN_DEFAULT_ALIGN_BYTES> MapVectorXs;
    //    MapVectorXs u = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, model_ref.nv, 1));
    //    {
    //      auto & u = internal_data.u;
    //      u.setZero();
    //      Eigen::Index row_id = 0;
    //      for (size_t constraint_id = 0; constraint_id < constraint_models_ref.size();
    //      ++constraint_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[constraint_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[constraint_id]);
    //        const auto csize = cmodel.size();
    //        const auto rhs_rows = rhs.middleRows(row_id, csize);
    //
    //        cmodel.jacobianTransposeMatrixProduct(model_ref, data_ref, cdata, rhs_rows, u,
    //        AddTo());
    //
    //        row_id += csize;
    //      }
    //    }

    // Backward sweep: propagate joint force contributions
    {
      //      for (auto & f : m_internal_data.f)
      //        f.setZero();
      //      auto & u = internal_data.u;
      //      u.setZero();

      typedef DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, data_ref, internal_data);
      for (JointIndex i = JointIndex(model_ref.njoints - 1); i > 0; --i)
      {
        Pass1::run(model_ref.joints[i], data_ref.joints[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass2;
      for (auto & motion : internal_data.a)
        motion.setZero();
      typename Pass2::ArgsType args2(model_ref, data_ref, internal_data);
      for (JointIndex i = 1; i < JointIndex(model_ref.njoints); ++i)
      {
        Pass2::run(model_ref.joints[i], data_ref.joints[i], args2);
      }
    }

    // Make a pass over the whole set of constraints to project back the accelerations onto the
    // joint
    mapJointSpaceToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, internal_data.a,
      internal_data.ddq, res, LocalFrameTag());

    // TODO(jcarpent): extend the code to operator on matrices
    //    {
    //      const auto & ddq = internal_data.ddq;
    //      Eigen::Index row_id = 0;
    //      for (size_t constraint_id = 0; constraint_id < constraint_models_ref.size();
    //      ++constraint_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[constraint_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[constraint_id]);
    //        const auto csize = cmodel.size();
    //
    //        cmodel.jacobianMatrixProduct(
    //          model_ref, data_ref, cdata, ddq, res.middleRows(row_id, csize));
    //
    //        row_id += csize;
    //      }
    //    }

    // Add damping contribution
    if (with_damping)
    {
      // res.array() += m_sum_compliance_damping.array() * rhs.array();
      m_sum_compliance_damping.template applyOnTheRight<::pinocchio::internal::add_assign_op>(
        rhs, res);
    }
    else
    {
      // take only compliance into account
      res.array() += m_compliance.array() * rhs.array();
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  template<typename MatrixLike>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat_) const
  {
    MatrixLike & mat = mat_.const_cast_derived();
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      mat.rows(), size(), "The input matrix does not match the size of the Delassus.");

    PINOCCHIO_THROW_IF(
      m_solve_in_place_dirty, std::logic_error,
      "The DelassusOperator has dirty quantities. Please call updateDecomposition() first.");

    const Model & model_ref = model();
    const Data & data_ref = data();
    const ConstraintModelVector & constraint_models_ref = constraint_models();
    const ConstraintDataVector & constraint_datas_ref = constraint_datas();
    auto & internal_data = this->m_internal_data;

    typedef Eigen::Map<VectorXs, EIGEN_DEFAULT_ALIGN_BYTES> MapVectorXs;
    MapVectorXs mat_tmp = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, size(), 1));

    // mat.array() *= m_sum_compliance_damping_inverse.array();
    m_sum_compliance_damping_inverse.template applyOnTheRight<::pinocchio::internal::assign_op>(
      mat, mat_tmp);
    mat = mat_tmp;

    // Make a pass over the whole set of constraints to add the contributions of constraint

    typedef Eigen::Map<VectorXs, EIGEN_DEFAULT_ALIGN_BYTES> MapVectorXs;
    MapVectorXs u = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, model_ref.nv, 1));
    // u and internal_data.of_augmented are reset by mapConstraintForcesToJointSpace
    mapConstraintForcesToJointSpace(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, mat,
      internal_data.of_augmented, u, WorldFrameTag());

    //    {
    //      u.setZero();
    //      Eigen::Index row_id = 0;
    //      for (size_t constraint_id = 0; constraint_id < constraint_models_ref.size();
    //      ++constraint_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[constraint_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[constraint_id]);
    //        const auto csize = cmodel.size();
    //        const auto mat_rows = mat.middleRows(row_id, csize);
    //
    //        cmodel.jacobianTransposeMatrixProduct(model_ref, data_ref, cdata, mat_rows, u,
    //        AddTo());
    //
    //        row_id += csize;
    //      }
    //    }

    const auto & augmented_mass_matrix_operator = this->getAugmentedMassMatrixOperator();
    augmented_mass_matrix_operator.solveInPlace(u, false);

    //    {
    //      Eigen::Index row_id = 0;
    //      for (size_t constraint_id = 0; constraint_id < constraint_models_ref.size();
    //      ++constraint_id)
    //      {
    //        const auto & cmodel =
    //          helper::get_ref(constraint_models_ref[constraint_id]);
    //        const auto & cdata =
    //          helper::get_ref(constraint_datas_ref[constraint_id]);
    //        const auto csize = cmodel.size();
    //
    //        cmodel.jacobianMatrixProduct(
    //          model_ref, data_ref, cdata, u, mat_tmp.middleRows(row_id, csize));
    //
    //        row_id += csize;
    //      }
    //    }

    // Make a pass over the whole set of constraints to project back the joint accelerations onto
    // the constraints
    mapJointSpaceToConstraintMotions(
      model_ref, data_ref, constraint_models_ref, constraint_datas_ref, internal_data.oa_augmented,
      u, mat_tmp, WorldFrameTag());

    // mat.noalias() -= m_sum_compliance_damping_inverse.asDiagonal() * mat_tmp;
    // m_sum_compliance_damping_inverse.template
    // applyOnTheRight<::pinocchio::internal::sub_assign_op>(mat_tmp,mat.noalias());
    m_sum_compliance_damping_inverse.template applyOnTheRight<::pinocchio::internal::sub_assign_op>(
      mat_tmp, mat); // TODO(jcarpent): fix me with proper noalias handling in applyOnTheRight
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder>
  template<typename MatrixLike>
  void DelassusOperatorRigidBodySystemsTpl<
    Scalar,
    Options,
    JointCollectionTpl,
    ConstraintModel,
    Holder>::AugmentedMassMatrixOperator::
    solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat_, bool reset_joint_force_vector) const
  {
    MatrixLike & mat = mat_.const_cast_derived();
    const auto & model_ref = m_self.model();
    const auto & data_ref = m_self.data();
    auto & internal_data =
      const_cast<DelassusOperatorRigidBodySystemsTpl &>(m_self).getInternalData();
    const auto & joint_elimination_order = data_ref.joint_elimination_order;

    if (reset_joint_force_vector)
    {
      for (auto & of_augmented : internal_data.of_augmented)
        of_augmented.setZero();
    }

    // Backward sweep: propagate joint force contributions
    {
      internal_data.u = mat;
      typedef AugmentedMassMatrixOperatorSolveInPlaceBackwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass1;
      typename Pass1::ArgsType args1(model_ref, data_ref, internal_data);
      for (const JointIndex i : joint_elimination_order)
      {
        Pass1::run(model_ref.joints[i], data_ref.joints_augmented[i], args1);
      }
    }

    // Forward sweep: compute joint accelerations
    {
      typedef AugmentedMassMatrixOperatorSolveInPlaceForwardPass<
        DelassusOperatorRigidBodySystemsTpl>
        Pass2;
      internal_data.oa_augmented[0].setZero();
      typename Pass2::ArgsType args2(model_ref, data_ref, internal_data);
      for (int it = int(joint_elimination_order.size()) - 1; it >= 0; it--)
      {
        const JointIndex i = joint_elimination_order[size_t(it)];
        Pass2::run(model_ref.joints[i], data_ref.joints_augmented[i], args2);
      }
    }

    mat = internal_data.ddq;
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hxx__
