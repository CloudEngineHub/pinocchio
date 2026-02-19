//
// Copyright (c) 2021-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_utils_hxx__
#define __pinocchio_algorithm_constraints_utils_hxx__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/utils/reference.hpp"

#include "pinocchio/math/matrix-block-type.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"

namespace pinocchio
{
  template<
    typename Scalar,
    int Options,
    int ForceOptions,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename ForceMatrix,
    class ForceAllocator,
    ReferenceFrame rf>
  void mapConstraintForcesToJointForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<ForceMatrix> & constraint_forces,
    std::vector<ForceTpl<Scalar, ForceOptions>, ForceAllocator> & joint_forces,
    ReferenceFrameTag<rf> reference_frame)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));

    const Eigen::Index constraint_size = getTotalConstraintResidualSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), constraint_size);

    for (auto & force : joint_forces)
      force.setZero();

    Eigen::Index row_id = 0;
    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
      const auto constraint_size = cmodel.residualSize();

      const auto constraint_force = constraint_forces.segment(row_id, constraint_size);
      cmodel.mapConstraintForceToJointForces(
        model, data, cdata, constraint_force, joint_forces, reference_frame);

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    int ForceOptions,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename ForceMatrix,
    class ForceAllocator,
    typename GeneralizedTorqueVector,
    ReferenceFrame rf>
  void mapConstraintForcesToJointSpace(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<ForceMatrix> & constraint_forces,
    std::vector<ForceTpl<Scalar, ForceOptions>, ForceAllocator> & joint_forces,
    const Eigen::MatrixBase<GeneralizedTorqueVector> & joint_torques_,
    ReferenceFrameTag<rf> reference_frame)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_torques_.size(), model.nv);

    const Eigen::Index constraint_size = getTotalConstraintResidualSize(constraint_models);
    assert(constraint_forces.rows() == constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), constraint_size);

    auto & joint_torques = joint_torques_.const_cast_derived();

    // Reset quantities
    joint_torques.setZero();
    for (auto & force : joint_forces)
      force.setZero();

    Eigen::Index row_id = 0;
    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
      const auto constraint_size = cmodel.residualSize();

      const auto constraint_force = constraint_forces.segment(row_id, constraint_size);
      cmodel.mapConstraintForceToJointSpace(
        model, data, cdata, constraint_force, joint_forces, joint_torques, reference_frame);

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    int MotionOptions,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    class MotionAllocator,
    typename MotionConstraintMatrix,
    ReferenceFrame rf>
  void mapJointMotionsToConstraintMotions(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const std::vector<MotionTpl<Scalar, MotionOptions>, MotionAllocator> & joint_motions,
    const Eigen::MatrixBase<MotionConstraintMatrix> & constraint_motions_,
    ReferenceFrameTag<rf> reference_frame)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.size(), size_t(model.njoints));

    auto & constraint_motions = constraint_motions_.const_cast_derived();
    const Eigen::Index constraint_size = getTotalConstraintResidualSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions.rows(), constraint_size);

    Eigen::Index row_id = 0;
    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
      const auto constraint_size = cmodel.residualSize();

      auto constraint_motion = constraint_motions.segment(row_id, constraint_size);
      cmodel.mapJointMotionsToConstraintMotion(
        model, data, cdata, joint_motions, constraint_motion, reference_frame);

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    int MotionOptions,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    class MotionAllocator,
    typename GeneralizedVelocityVector,
    typename MotionConstraintMatrix,
    ReferenceFrame rf>
  void mapJointSpaceToConstraintMotions(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const std::vector<MotionTpl<Scalar, MotionOptions>, MotionAllocator> & joint_motions,
    const Eigen::MatrixBase<GeneralizedVelocityVector> & generalized_velocity,
    const Eigen::MatrixBase<MotionConstraintMatrix> & constraint_motions_,
    ReferenceFrameTag<rf> reference_frame)
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_models.size(), constraint_datas.size());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.size(), size_t(model.njoints));
    PINOCCHIO_CHECK_ARGUMENT_SIZE(generalized_velocity.size(), model.nv);

    auto & constraint_motions = constraint_motions_.const_cast_derived();
    const Eigen::Index total_constraint_size = getTotalConstraintResidualSize(constraint_models);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions.rows(), total_constraint_size);

    Eigen::Index row_id = 0;
    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
      const auto constraint_size = cmodel.residualSize();

      auto constraint_motion = constraint_motions.segment(row_id, constraint_size);
      cmodel.mapJointSpaceToConstraintMotion(
        model, data, cdata, joint_motions, generalized_velocity, constraint_motion,
        reference_frame);

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModel,
    typename ConstraintData,
    typename JacobianMatrixLike>
  void getConstraintJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintModelBase<ConstraintModel> & constraint_model_,
    const ConstraintDataBase<ConstraintData> & constraint_data_,
    const Eigen::MatrixBase<JacobianMatrixLike> & J_)
  {
    JacobianMatrixLike & J = J_.const_cast_derived();
    const auto & constraint_model = helper::get_ref(constraint_model_.derived());
    const auto & constraint_data = helper::get_ref(constraint_data_.derived());

    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_model.residualSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    constraint_model.jacobian(model, data, constraint_data, J);
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename DynamicMatrixLike>
  void getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<DynamicMatrixLike> & J_)
  {
    const Eigen::Index constraint_size = getTotalConstraintResidualSize(constraint_models);
    assert(J_.rows() == constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(J_.cols(), model.nv);

    assert(model.check(data) && "data is not consistent with model.");
    assert(model.check(MimicChecker()) && "Function does not support mimic joints");

    auto & J = J_.const_cast_derived();
    Eigen::Index row_id = 0;
    for (size_t k = 0; k < constraint_models.size(); ++k)
    {
      const auto & cmodel = helper::get_ref(constraint_models[k]);
      const auto & cdata = helper::get_ref(constraint_datas[k]);

      const auto csize = cmodel.residualSize();
      getConstraintJacobian(model, data, cmodel, cdata, J.middleRows(row_id, csize));

      row_id += csize;
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  typename DataTpl<Scalar, Options, JointCollectionTpl>::MatrixXs getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
    typedef typename Data::MatrixXs ReturnType;

    const auto constraint_size = getTotalConstraintResidualSize(constraint_models);

    ReturnType res = ReturnType::Zero(constraint_size, model.nv);
    getConstraintsJacobian(model, data, constraint_models, constraint_datas, res);

    return res;
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType,
    AssignmentOperatorType op>
  void evalConstraintJacobianMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res_,
    AssignmentOperatorTag<op> aot)
  {
    PINOCCHIO_UNUSED_VARIABLE(aot);
    const Eigen::Index constraint_size = getTotalConstraintResidualSize(constraint_models);
    auto & res = res_.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(rhs.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.cols(), rhs.cols());

    using aot_internal = std::conditional_t<
      std::is_same<AssignmentOperatorTag<op>, SetTo>::value, AddTo, AssignmentOperatorTag<op>>;
    if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
    {
      res.setZero();
    }

    Eigen::Index row_id = 0;
    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
      const auto constraint_size = cmodel.residualSize();

      auto res_block = res.middleRows(row_id, constraint_size);
      cmodel.jacobianMatrixProduct(model, data, cdata, rhs, res_block, aot_internal());

      row_id += constraint_size;
    }
  }

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType,
    AssignmentOperatorType op>
  void evalConstraintJacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res_,
    AssignmentOperatorTag<op> aot)
  {
    PINOCCHIO_UNUSED_VARIABLE(aot);
    const Eigen::Index constraint_size = getTotalConstraintResidualSize(constraint_models);
    ResultMatrixType & res = res_.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(rhs.rows(), constraint_size);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res_.cols(), rhs.cols());

    Eigen::Index row_id = 0;

    using aot_internal = std::conditional_t<
      std::is_same<AssignmentOperatorTag<op>, SetTo>::value, AddTo, AssignmentOperatorTag<op>>;
    if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
    {
      res.setZero();
    }

    for (size_t constraint_id = 0; constraint_id < constraint_models.size(); ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
      const auto constraint_size = cmodel.residualSize();

      const auto rhs_block = rhs.middleRows(row_id, constraint_size);
      cmodel.jacobianTransposeMatrixProduct(model, data, cdata, rhs_block, res, aot_internal());

      row_id += constraint_size;
    }
  }

  namespace internal
  {
    /**
     * @brief ComputeBlockDiagonalPatternImpl functor
     * Computes the block pattern for a specific constraint.
     */

    // Default implementation of computing the block diagonal matrix block_infos
    // for this constraint.
    template<typename ConstraintModel, class = void>
    struct ComputeBlockDiagonalPatternImpl
    {
      template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
      static void run(
        const ConstraintModel & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        PINOCCHIO_UNUSED_VARIABLE(cmodel);
        PINOCCHIO_UNUSED_VARIABLE(block_infos);
        PINOCCHIO_UNUSED_VARIABLE(dispatcher);

        static_assert(
          false, "ComputeBlockDiagonalPatternImpl not implemented for this constraint.");
      }
    };

    // Specialization of compute block_infos for FrameAnchor.
    template<typename Scalar, int Options>
    struct ComputeBlockDiagonalPatternImpl<FrameAnchorConstraintModelTpl<Scalar, Options>>
    {
      typedef FrameAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;

      template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
      static void run(
        const ConstraintModel & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        PINOCCHIO_UNUSED_VARIABLE(cmodel);
        PINOCCHIO_UNUSED_VARIABLE(dispatcher);
        assert(ConstraintModel::SymmetricConeSize != Eigen::Dynamic);

        // equality constraint -> prox term -> Scalar Identity
        block_infos.emplace_back(MatrixBlockType::ScalarIdentity, ConstraintModel::Size);
      }
    };

    // Specialization of compute block_infos for PointAnchor.
    template<typename Scalar, int Options>
    struct ComputeBlockDiagonalPatternImpl<PointAnchorConstraintModelTpl<Scalar, Options>>
    {
      typedef PointAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;

      template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
      static void run(
        const ConstraintModel & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        PINOCCHIO_UNUSED_VARIABLE(cmodel);
        PINOCCHIO_UNUSED_VARIABLE(dispatcher);
        assert(ConstraintModel::SymmetricConeSize != Eigen::Dynamic);

        // equality constraint -> prox term -> Scalar Identity
        block_infos.emplace_back(MatrixBlockType::ScalarIdentity, ConstraintModel::Size);
      }
    };

    // Specialization of compute block_infos for PointContact.
    template<typename Scalar, int Options>
    struct ComputeBlockDiagonalPatternImpl<PointContactConstraintModelTpl<Scalar, Options>>
    {
      typedef PointContactConstraintModelTpl<Scalar, Options> ConstraintModel;

      template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
      static void run(
        const ConstraintModel & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        PINOCCHIO_UNUSED_VARIABLE(cmodel);
        PINOCCHIO_UNUSED_VARIABLE(dispatcher);
        assert(ConstraintModel::SymmetricConeSize != Eigen::Dynamic);

        if constexpr (std::is_same<BlockDiagonalDispatcherTag<op>, DiagonalDispatcher>::value)
        {
          block_infos.emplace_back(MatrixBlockType::ScalarIdentity, ConstraintModel::Size);
        }

        if constexpr (std::is_same<
                        BlockDiagonalDispatcherTag<op>, IPMBlockDiagonalDispatcher>::value)
        {
          // for ipm, point contact -> 3x3 block
          block_infos.emplace_back(MatrixBlockType::Plain, ConstraintModel::Size);
        }
      }
    };

    // Specialization of compute block_infos for JointLimit.
    template<typename Scalar, int Options>
    struct ComputeBlockDiagonalPatternImpl<JointLimitConstraintModelTpl<Scalar, Options>>
    {
      typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;

      template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
      static void run(
        const ConstraintModel & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        PINOCCHIO_UNUSED_VARIABLE(dispatcher);

        if constexpr (std::is_same<BlockDiagonalDispatcherTag<op>, DiagonalDispatcher>::value)
        {
          block_infos.emplace_back(MatrixBlockType::ScalarIdentity, cmodel.residualSize());
        }

        if constexpr (std::is_same<
                        BlockDiagonalDispatcherTag<op>, IPMBlockDiagonalDispatcher>::value)
        {
          block_infos.emplace_back(MatrixBlockType::Diagonal, cmodel.residualSize());
        }
      }
    };

    // Specialization of compute block_infos for JointFriction.
    template<typename Scalar, int Options>
    struct ComputeBlockDiagonalPatternImpl<JointFrictionConstraintModelTpl<Scalar, Options>>
    {
      typedef JointFrictionConstraintModelTpl<Scalar, Options> ConstraintModel;

      template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
      static void run(
        const ConstraintModel & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        PINOCCHIO_UNUSED_VARIABLE(dispatcher);

        if constexpr (std::is_same<BlockDiagonalDispatcherTag<op>, DiagonalDispatcher>::value)
        {
          block_infos.emplace_back(MatrixBlockType::ScalarIdentity, cmodel.residualSize());
        }

        if constexpr (std::is_same<
                        BlockDiagonalDispatcherTag<op>, IPMBlockDiagonalDispatcher>::value)
        {
          block_infos.emplace_back(MatrixBlockType::Diagonal, cmodel.residualSize());
        }
      }
    };

    /**
     * @brief ComputeBlockDiagonalPatternVisitor visitor
     */
    template<typename BlockInfoVector, BlockDiagonalDispatcherType op>
    struct ComputeBlockDiagonalPatternVisitor
    : visitors::ConstraintUnaryVisitorBase<ComputeBlockDiagonalPatternVisitor<BlockInfoVector, op>>
    {
      typedef boost::fusion::vector<BlockInfoVector &, BlockDiagonalDispatcherTag<op>> ArgsType;

      typedef visitors::ConstraintUnaryVisitorBase<
        ComputeBlockDiagonalPatternVisitor<BlockInfoVector, op>>
        Base;
      using Base::run;

      template<typename ConstraintModel>
      static void algo(
        const ConstraintModelBase<ConstraintModel> & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        typedef ComputeBlockDiagonalPatternImpl<ConstraintModel> Impl;
        Impl::run(cmodel.derived(), block_infos, dispatcher);
      }

      template<typename ConstraintModel>
      static void run(
        const pinocchio::ConstraintModelBase<ConstraintModel> & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        algo(cmodel.derived(), block_infos, dispatcher);
      }

      template<
        typename Scalar,
        int Options,
        template<typename S, int O> class ConstraintCollectionTpl>
      static void run(
        const pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
        BlockInfoVector & block_infos,
        BlockDiagonalDispatcherTag<op> dispatcher)
      {
        ArgsType args(block_infos, dispatcher);
        run(cmodel.derived(), args);
      }
    }; // struct ComputeBlockDiagonalPatternVisitor
  } // namespace internal

  template<
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename BlockDiagonalElement,
    BlockDiagonalDispatcherType op>
  void computeBlockDiagonalPattern(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    std::vector<BlockDiagonalElement> & block_diagonal_infos,
    BlockDiagonalDispatcherTag<op> dispatcher)
  {
    block_diagonal_infos.clear();
    block_diagonal_infos.reserve(constraint_models.size());

    for (std::size_t i = 0; i < constraint_models.size(); ++i)
    {
      const auto & cmodel = helper::get_ref(constraint_models[i]);

      typedef internal::ComputeBlockDiagonalPatternVisitor<decltype(block_diagonal_infos), op> Algo;
      Algo::run(cmodel, block_diagonal_infos, dispatcher);
    }
  }

  template<
    typename ConstraintModel,
    typename ConstraintModelAllocator,
    typename Scalar,
    int Options,
    std::size_t Alignment>
  void constructPositiveDefiniteBlockDiagonalMatrix(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    BlockDiagonalMatrixTpl<Scalar, Options, Alignment> & block_diagonal_matrix)
  {
    typedef BlockDiagonalMatrixTpl<Scalar, Options, Alignment> BlockDiagonalMatrix;
    typedef typename BlockDiagonalMatrix::MatrixBlockElement MatrixBlockElement;
    std::vector<MatrixBlockElement> block_diagonal_infos;
    computeBlockDiagonalPattern(
      constraint_models, block_diagonal_infos, IPMBlockDiagonalDispatcher());
    block_diagonal_matrix.rebuild(block_diagonal_infos);
    for (auto & block : block_diagonal_matrix.blocks())
    {
      block.setRandomPD();
    }
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_utils_hxx__
