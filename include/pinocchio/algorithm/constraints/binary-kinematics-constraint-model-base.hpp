//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_relative_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraints_relative_constraint_model_base_hpp__

#include <algorithm>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/kinematics-constraint-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct BinaryKinematicsConstraintBase
  : KinematicsConstraintModelBase<Derived>
  , ConstraintModelCommonParameters<Derived>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };

    typedef KinematicsConstraintModelBase<Derived> KinematicsBase;
    typedef ConstraintModelCommonParameters<Derived> BaseCommonParameters;
    typedef ConstraintModelBase<Derived> RootBase;

    template<typename OtherDerived>
    friend struct BinaryKinematicsConstraintBase;

    using KinematicsBase::joint1_id;
    using KinematicsBase::joint2_id;
    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef typename traits<Derived>::VectorConstraintSize VectorConstraintSize;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    using RootBase::maxResidualSize;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // CRTP related ------------------

    /// \brief Cast to Base
    KinematicsBase & base()
    {
      return static_cast<KinematicsBase &>(*this);
    }

    /// \brief Const cast to Base
    const KinematicsBase & base() const
    {
      return static_cast<const KinematicsBase &>(*this);
    }

    /// \brief Cast to BaseCommonParameters
    BaseCommonParameters & base_common_parameters()
    {
      return static_cast<BaseCommonParameters &>(*this);
    }

    /// \brief Const cast to BaseCommonParameters
    const BaseCommonParameters & base_common_parameters() const
    {
      return static_cast<const BaseCommonParameters &>(*this);
    }

    // Constructors ------------------

    /// \brief Default constructor
    BinaryKinematicsConstraintBase()
    : joint1_placement(SE3::Identity())
    , joint2_placement(SE3::Identity())
    , desired_constraint_offset(VectorConstraintSize::Zero())
    , desired_constraint_velocity(VectorConstraintSize::Zero())
    , desired_constraint_acceleration(VectorConstraintSize::Zero())
    , nv(-1)
    , depth_joint1(0)
    , depth_joint2(0)
    {
    }

    /// \brief Full constructor
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : KinematicsBase(model, joint1_id, joint2_id)
    , joint1_placement(joint1_placement)
    , joint2_placement(joint2_placement)
    , desired_constraint_offset(VectorConstraintSize::Zero())
    , desired_constraint_velocity(VectorConstraintSize::Zero())
    , desired_constraint_acceleration(VectorConstraintSize::Zero())
    , colwise_joint1_sparsity(model.nv)
    , colwise_joint2_sparsity(model.nv)
    , nv(-1)
    , depth_joint1(0)
    , depth_joint2(0)
    {
      init(model);
    }

    /// \brief Constructor with only joint1, relative placement is identity.
    /// joint2 defaults to 0 with identity relative placement.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : BinaryKinematicsConstraintBase(model, joint1_id, SE3::Identity(), 0, SE3::Identity())
    {
    }

    /// \brief Constructor with only joint1 and relative placement to joint1.
    /// joint2 defaults to 0 with identity relative placement.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : BinaryKinematicsConstraintBase(model, joint1_id, joint1_placement, 0, SE3::Identity())
    {
    }

    /// \brief Constructor with only joint1 and joint2 ids.
    /// Relative placements are identity.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : BinaryKinematicsConstraintBase(model, joint1_id, SE3::Identity(), joint2_id, SE3::Identity())
    {
    }

    // Operators ---------------------

    /// \brief Comparison operator.
    template<typename OtherDerived>
    bool operator==(const BinaryKinematicsConstraintBase<OtherDerived> & other) const
    {
      if (this == &other)
        return true;

      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && joint1_id == other.joint1_id && joint2_id == other.joint2_id
             && joint1_placement == other.joint1_placement
             && joint2_placement == other.joint2_placement && nv == other.nv
             && desired_constraint_offset == other.desired_constraint_offset
             && desired_constraint_velocity == other.desired_constraint_velocity
             && desired_constraint_acceleration == other.desired_constraint_acceleration
             && colwise_joint1_sparsity == other.colwise_joint1_sparsity
             && colwise_joint2_sparsity == other.colwise_joint2_sparsity
             && joint1_span_indexes == other.joint1_span_indexes
             && joint2_span_indexes == other.joint2_span_indexes
             && depth_joint1 == other.depth_joint1 && depth_joint2 == other.depth_joint2
             && colwise_sparsity == other.colwise_sparsity
             && colwise_span_indexes == other.colwise_span_indexes;
    }

    /// \brief Comparison operator.
    template<typename OtherDerived>
    bool operator!=(const BinaryKinematicsConstraintBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

    /// \brief Cast to NewScalar.
    template<typename NewScalar, typename OtherDerived>
    void cast(BinaryKinematicsConstraintBase<OtherDerived> & res) const
    {
      KinematicsBase::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.joint1_id = joint1_id;
      res.joint2_id = joint2_id;
      res.joint1_placement = joint1_placement.template cast<NewScalar>();
      res.joint2_placement = joint2_placement.template cast<NewScalar>();
      res.desired_constraint_offset = desired_constraint_offset.template cast<NewScalar>();
      res.desired_constraint_velocity = desired_constraint_velocity.template cast<NewScalar>();
      res.desired_constraint_acceleration =
        desired_constraint_acceleration.template cast<NewScalar>();
      res.colwise_joint1_sparsity = colwise_joint1_sparsity;
      res.colwise_joint2_sparsity = colwise_joint2_sparsity;
      res.joint1_span_indexes = joint1_span_indexes;
      res.joint2_span_indexes = joint2_span_indexes;
      res.colwise_sparsity = colwise_sparsity;
      res.colwise_span_indexes = colwise_span_indexes;
      res.nv = nv;
      res.depth_joint1 = depth_joint1;
      res.depth_joint2 = depth_joint2;
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc Base::getRowSparsityPattern
    template<template<typename, int> class JointCollectionTpl>
    const BooleanVector & getRowSparsityPatternImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < maxResidualSize());
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return colwise_sparsity;
    }

    /// \copydoc Base::getRowIndexes
    template<template<typename, int> class JointCollectionTpl>
    const EigenIndexVector & getRowIndexesImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < maxResidualSize());
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return colwise_span_indexes;
    }

  protected:
    // ------------------------------
    // PROTECTED METHODS
    // ------------------------------

    /// \brief Initialize the constraint model based on the kinematics loop of model.
    template<template<typename, int> class JointCollectionTpl>
    void init(const ModelTpl<Scalar, Options, JointCollectionTpl> & model);

  public:
    // ------------------------------
    // MEMBERS
    // ------------------------------

    /// \brief Position of attached point with respect to the frame of joint1.
    SE3 joint1_placement;

    /// \brief Position of attached point with respect to the frame of joint2.
    SE3 joint2_placement;

    /// \brief Desired constraint shift at position level.
    VectorConstraintSize desired_constraint_offset;

    /// \brief Desired constraint velocity at velocity level.
    VectorConstraintSize desired_constraint_velocity;

    /// \brief Desired constraint acceleration at acceleration level.
    VectorConstraintSize desired_constraint_acceleration;

    /// \brief Column-wise sparsity pattern associated with joint 1.
    BooleanVector colwise_joint1_sparsity;

    /// \brief Column-wise sparsity pattern associated with joint 2.
    BooleanVector colwise_joint2_sparsity;

    /// \brief Joint-wise span indexes associated with joint 1.
    EigenIndexVector joint1_span_indexes;

    /// \brief Joint-wise span indexes associated with joint 2.
    EigenIndexVector joint2_span_indexes;

    /// \brief Sparsity pattern associated to the constraint.
    BooleanVector colwise_sparsity;

    /// \brief Indexes of the columns spanned by the constraints.
    EigenIndexVector colwise_span_indexes;

    /// \brief Dimensions of the model.
    int nv;

    /// \brief Depth of the kinematic tree for joint1 and joint2.
    size_t depth_joint1, depth_joint2;

  protected:
    using BaseCommonParameters::m_baumgarte_parameters;
    using BaseCommonParameters::m_compliance;
  }; // struct BinaryKinematicsConstraintBase<Derived>

} // namespace pinocchio

#include "pinocchio/algorithm/constraints/binary-kinematics-constraint-model-base.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_relative_constraint_model_base_hpp__
