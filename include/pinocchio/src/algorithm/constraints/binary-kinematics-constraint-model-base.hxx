//
// Copyright (c) 2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/algorithm/constraints.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/algorithm/constraints.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  // --------------------------------------------------------------
  // Declaration
  // --------------------------------------------------------------
  template<typename Derived>
  struct BinaryKinematicsConstraintModelBase;

  // --------------------------------------------------------------
  // Helpers
  // --------------------------------------------------------------
  template<typename Derived>
  using enable_if_binarykinematic_model_t =
    std::enable_if_t<std::is_base_of_v<BinaryKinematicsConstraintModelBase<Derived>, Derived>>;

  // --------------------------------------------------------------
  // Struct
  // --------------------------------------------------------------
  template<typename Derived>
  struct BinaryKinematicsConstraintModelBase
  : KinematicsConstraintModelBase<Derived>
  , ConstraintModelCommonParameters<Derived>
  {
    // --------------------------------------------------------------
    // Type defs
    // --------------------------------------------------------------
    // CRTP related types -------------------------------------------
    typedef KinematicsConstraintModelBase<Derived> Base;
    typedef ConstraintModelCommonParameters<Derived> BaseCommonParameters;
    typedef ConstraintModelBase<Derived> RootBase;

    // Retrieving traits --------------------------------------------
    typedef typename traits<Derived>::ConstraintModel ConstraintModel;
    typedef typename traits<Derived>::ConstraintData ConstraintData;

    typedef typename traits<Derived>::Scalar Scalar;
    static constexpr int Options = traits<Derived>::Options;

    static constexpr int Size = traits<Derived>::Size;

    typedef typename traits<Derived>::ResidualVectorType ResidualVectorType;
    typedef typename traits<Derived>::JacobianMatrixType JacobianMatrixType;
    typedef typename traits<Derived>::BaumgarteCorrectorParameters BaumgarteCorrectorParameters;

    // Friendship ---------------------------------------------------
    template<typename OtherDerived>
    friend struct BinaryKinematicsConstraintModelBase;

    // Base usage --------------------------------------------------
    using RootBase::derived;
    using RootBase::residualSize;
    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    // Useful types ------------------------------------------------
    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef ForceTpl<Scalar, Options> Force;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef Eigen::Matrix<Scalar, 3, 6, Options> Matrix36;
    typedef Eigen::Matrix<Scalar, Size, 6, Options> MatrixSize6;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // CRTP related ------------------

    /// \brief Cast to Base
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    /// \brief Const cast to Base
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
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

  protected:
    /// \brief Default constructor
    BinaryKinematicsConstraintModelBase()
    : joint1_id(0)
    , joint2_id(0)
    , joint1_placement(SE3::Identity())
    , joint2_placement(SE3::Identity())
    , desired_constraint_offset(ResidualVectorType::Zero(residualSize()))
    , desired_constraint_velocity(ResidualVectorType::Zero(residualSize()))
    , desired_constraint_acceleration(ResidualVectorType::Zero(residualSize()))
    , nv(-1)
    , depth_joint1(0)
    , depth_joint2(0)
    {
    }

    /// \brief Constructor with only model.
    /// Relative placements are identity and joints are 0.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    : BinaryKinematicsConstraintModelBase(model, 0, SE3::Identity(), 0, SE3::Identity())
    {
    }

    /// \brief Full constructor
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : Base(model)
    , joint1_id(joint1_id)
    , joint2_id(joint2_id)
    , joint1_placement(joint1_placement)
    , joint2_placement(joint2_placement)
    , desired_constraint_offset(ResidualVectorType::Zero(residualSize()))
    , desired_constraint_velocity(ResidualVectorType::Zero(residualSize()))
    , desired_constraint_acceleration(ResidualVectorType::Zero(residualSize()))
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
    BinaryKinematicsConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : BinaryKinematicsConstraintModelBase(model, joint1_id, SE3::Identity(), 0, SE3::Identity())
    {
    }

    /// \brief Constructor with only joint1 and relative placement to joint1.
    /// joint2 defaults to 0 with identity relative placement.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : BinaryKinematicsConstraintModelBase(model, joint1_id, joint1_placement, 0, SE3::Identity())
    {
    }

    /// \brief Constructor with only joint1 and joint2 ids.
    /// Relative placements are identity.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    BinaryKinematicsConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : BinaryKinematicsConstraintModelBase(
        model, joint1_id, SE3::Identity(), joint2_id, SE3::Identity())
    {
    }

    // Operators ---------------------

  public:
    /// \brief Cast to NewScalar.
    template<typename NewScalar, typename OtherDerived>
    void cast(BinaryKinematicsConstraintModelBase<OtherDerived> & res) const
    {
      Base::cast(res);
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

    /// \brief Comparison operator.
    template<typename OtherDerived>
    bool operator==(const BinaryKinematicsConstraintModelBase<OtherDerived> & other) const
    {
      if (this == &other)
        return true;

      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && joint1_id == other.joint1_id && joint2_id == other.joint2_id
             && joint1_placement == other.joint1_placement
             && joint2_placement == other.joint2_placement
             && desired_constraint_offset == other.desired_constraint_offset
             && desired_constraint_velocity == other.desired_constraint_velocity
             && desired_constraint_acceleration == other.desired_constraint_acceleration
             && colwise_joint1_sparsity == other.colwise_joint1_sparsity
             && colwise_joint2_sparsity == other.colwise_joint2_sparsity
             && joint1_span_indexes == other.joint1_span_indexes
             && joint2_span_indexes == other.joint2_span_indexes
             && colwise_sparsity == other.colwise_sparsity
             && colwise_span_indexes == other.colwise_span_indexes && nv == other.nv
             && depth_joint1 == other.depth_joint1 && depth_joint2 == other.depth_joint2;
    }

    /// \brief Comparison operator.
    template<typename OtherDerived>
    bool operator!=(const BinaryKinematicsConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

    // Binary related ----------------

    /// \brief Returns the constraint projector associated with joint 1.
    /// This matrix transforms a spatial velocity expressed at the origin to the first component of
    /// the constraint associated with joint 1.
    template<ReferenceFrame rf>
    MatrixSize6 getA1(const ConstraintData & cdata, ReferenceFrameTag<rf> rft) const
    {
      return derived().getA1Impl(cdata, rft);
    }

    /// \brief Returns the constraint projector associated with joint 2.
    /// This matrix transforms a spatial velocity expressed at the origin to the first component of
    /// the constraint associated with joint 2.
    template<ReferenceFrame rf>
    MatrixSize6 getA2(const ConstraintData & cdata, ReferenceFrameTag<rf> rft) const
    {
      return derived().getA2Impl(cdata, rft);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc Base::getRowSparsityPattern
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    const BooleanVector & getRowSparsityPatternImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < residualSize());
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return colwise_sparsity;
    }

    /// \copydoc Base::getRowIndexes
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    const EigenIndexVector & getRowIndexesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < residualSize());
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
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void init(const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model);

  public:
    // ------------------------------
    // MEMBERS
    // ------------------------------

    /// \brief Index of the first joint in the model tree
    JointIndex joint1_id;

    /// \brief Index of the second joint in the model tree
    JointIndex joint2_id;

    /// \brief Position of attached point with respect to the frame of joint1.
    SE3 joint1_placement;

    /// \brief Position of attached point with respect to the frame of joint2.
    SE3 joint2_placement;

    /// \brief Desired constraint shift at position level.
    ResidualVectorType desired_constraint_offset;

    /// \brief Desired constraint velocity at velocity level.
    ResidualVectorType desired_constraint_velocity;

    /// \brief Desired constraint acceleration at acceleration level.
    ResidualVectorType desired_constraint_acceleration;

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
  }; // struct BinaryKinematicsConstraintModelBase

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
