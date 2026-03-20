//
// Copyright (c) 2024-2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/constraints.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/constraints.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  // --------------------------------------------------------------
  // Cast
  // --------------------------------------------------------------
  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, JointFrictionConstraintModelTpl<Scalar, Options>>
  {
    typedef JointFrictionConstraintModelTpl<NewScalar, Options> type;
  };

  // --------------------------------------------------------------
  // Traits
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct traits<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef JointFrictionConstraintModelTpl<_Scalar, _Options> ConstraintModel;
    typedef JointFrictionConstraintDataTpl<_Scalar, _Options> ConstraintData;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    // --------------------------------------------------------------
    // Traits characterizing the constraints
    // --------------------------------------------------------------
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    static constexpr ConstraintSizeType constraint_size_type = ConstraintSizeType::CONSTANT;

    static constexpr bool has_baumgarte_corrector = false;
    static constexpr bool has_set = true;
    static constexpr bool is_inequality_constraint = true;

    // --------------------------------------------------------------
    // Traits for associated struct and sizes
    // --------------------------------------------------------------
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;
    typedef NonNegativeOrthantJordanOperationTpl<Scalar, Options> JordanOperation;
    typedef boost::blank BaumgarteCorrectorParameters;

    static constexpr int Size = Eigen::Dynamic;
    static constexpr int SymmetricConeSize = JordanOperation::ConeSize;
    static constexpr int SymmetricConeScalingSize = JordanOperation::ConeScalingSize;

    // --------------------------------------------------------------
    // Traits that are helper for Eigen types
    // --------------------------------------------------------------
    typedef Eigen::Matrix<Scalar, Size, 1, Options> ResidualVectorType;
    typedef Eigen::Matrix<Scalar, Size, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Eigen::Matrix<Scalar, SymmetricConeSize, 1, Options> ConeVectorType;
    typedef Eigen::Matrix<Scalar, SymmetricConeScalingSize, 1, Options> ConeScalingVectorType;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXs;

    // Template to generate type
    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::
        Matrix<Scalar, Size, InputMatrixPlain::ColsAtCompileTime, InputMatrixPlain::Options>
          type;
    };

    template<typename InputMatrix>
    struct JacobianTransposeMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Eigen::Dynamic,
        InputMatrixPlain::ColsAtCompileTime,
        InputMatrixPlain::Options>
        type;
    };
  };

  template<typename _Scalar, int _Options>
  struct traits<JointFrictionConstraintDataTpl<_Scalar, _Options>>
  : traits<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  {
  };

  // --------------------------------------------------------------
  // Struct
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct JointFrictionConstraintModelTpl
  : JointWiseConstraintModelBase<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Type defs
    // --------------------------------------------------------------
    // CRTP related types -------------------------------------------
    typedef JointFrictionConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelCommonParameters<Self> BaseCommonParameters;
    typedef ConstraintModelBase<JointFrictionConstraintModelTpl> RootBase;

    // Retrieving traits --------------------------------------------
    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintData ConstraintData;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    static constexpr ConstraintSizeType constraint_size_type = traits<Self>::constraint_size_type;

    static constexpr bool has_baumgarte_corrector = traits<Self>::has_baumgarte_corrector;

    typedef typename traits<Self>::ConstraintSet ConstraintSet;
    typedef typename traits<Self>::JordanOperation JordanOperation;
    typedef typename traits<Self>::BaumgarteCorrectorParameters BaumgarteCorrectorParameters;

    static constexpr int Size = traits<Self>::Size;
    static constexpr int SymmetricConeSize = traits<Self>::SymmetricConeSize;
    static constexpr int SymmetricConeScalingSize = traits<Self>::SymmetricConeScalingSize;

    typedef typename traits<Self>::ResidualVectorType ResidualVectorType;
    typedef typename traits<Self>::JacobianMatrixType JacobianMatrixType;
    typedef typename traits<Self>::ConeVectorType ConeVectorType;
    typedef typename traits<Self>::ConeScalingVectorType ConeScalingVectorType;

    // Friendship ---------------------------------------------------
    template<typename NewScalar, int NewOptions>
    friend struct JointFrictionConstraintModelTpl;

    // Base usage ---------------------------------------------------
    using RootBase::classname;
    using RootBase::jacobianMatrixProduct;
    using RootBase::jacobianTransposeMatrixProduct;
    using RootBase::residualSize;
    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    // Useful types ------------------------------------------------
    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectorOfEigenIndexVector;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;

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

    /// \brief Default constructor
    JointFrictionConstraintModelTpl()
    {
    }

    /// \brief Constructor from model only.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    JointFrictionConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    {
      size_t n_joints = model.joints.size();
      m_active_joints.reserve(n_joints);
      for (size_t i = 0; i < n_joints; ++i)
      {
        m_active_joints.push_back(static_cast<JointIndex>(i));
      }
      init(model);
    }

    /// \brief Constructor from model and active_joints.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    JointFrictionConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndexVector & active_joints)
    : Base(model)
    , m_active_joints(active_joints)
    {
      init(model);
    }

    // Operators ---------------------

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, JointFrictionConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, JointFrictionConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.m_active_joints = m_active_joints;
      res.m_active_dofs = m_active_dofs;
      res.m_active_joint_ids = m_active_joint_ids;
      res.m_friction_lower_limit = m_friction_lower_limit.template cast<NewScalar>();
      res.m_friction_upper_limit = m_friction_upper_limit.template cast<NewScalar>();
      return res;
    }

    ///
    /// \brief Comparison operator
    ///
    /// \param[in] other Other JointFrictionConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributes
    /// must be the same).
    ///
    bool operator==(const JointFrictionConstraintModelTpl & other) const
    {
      if (this == &other)
        return true;
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && m_active_joints == other.m_active_joints && m_active_dofs == other.m_active_dofs
             && m_active_joint_ids == other.m_active_joint_ids
             && m_friction_lower_limit == other.m_friction_lower_limit
             && m_friction_upper_limit == other.m_friction_upper_limit;
    }

    /// \brief Comparison operator
    bool operator!=(const JointFrictionConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    // Managing methods --------------

    /// \brief Returns the vector of active joints
    const JointIndexVector & getActiveJoints() const
    {
      return m_active_joints;
    }

    /// \brief Returns the vector of active rows
    const EigenIndexVector & getActiveDofs() const
    {
      return m_active_dofs;
    }

    /// \brief Returns a const reference to `lower_friction_limit`
    /// \note The upper/lower friction limits units should be coherent with the algos.
    /// If an algo works with forces, limits should be expressed in Newtons.
    /// If an algo works with impulses, limits should be expressed in Newtons * Time.
    /// Typically, constraint solvers work on impulses.
    const VectorXs & getFrictionLowerLimit() const
    {
      return m_friction_lower_limit;
    }

    /// \brief Set the lower friction limit.
    /// \note The upper/lower friction limits units should be coherent with the algos.
    /// If an algo works with forces, limits should be expressed in Newtons.
    /// If an algo works with impulses, limits should be expressed in Newtons * Time.
    /// Typically, constraint solvers work on impulses.
    template<typename VectorLike>
    void setFrictionLowerLimit(const Eigen::MatrixBase<VectorLike> & lb)
    {
      PINOCCHIO_THROW_IF(
        lb.size() != residualSize(), std::runtime_error, "lb should be the same as size()");
      m_friction_lower_limit = lb;
    }

    /// \brief Returns a const reference to `upper_friction_limit`
    /// \note The upper/lower friction limits units should be coherent with the algos.
    /// If an algo works with forces, limits should be expressed in Newtons.
    /// If an algo works with impulses, limits should be expressed in Newtons * Time.
    /// Typically, constraint solvers work on impulses.
    const VectorXs & getFrictionUpperLimit() const
    {
      return m_friction_upper_limit;
    }

    /// \brief Set the upper friction limit.
    /// \note The upper/lower friction limits units should be coherent with the algos.
    /// If an algo works with forces, limits should be expressed in Newtons.
    /// If an algo works with impulses, limits should be expressed in Newtons * Time.
    /// Typically, constraint solvers work on impulses.
    template<typename VectorLike>
    void setFrictionUpperLimit(const Eigen::MatrixBase<VectorLike> & ub)
    {
      PINOCCHIO_THROW_IF(
        ub.size() != residualSize(), std::runtime_error, "ub should be the same as size()");
      m_friction_upper_limit = ub;
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    // General -----------------------

    /// \copydoc RootBase::classname
    static std::string classnameImpl()
    {
      return std::string("JointFrictionConstraintModel");
    }

    /// \copydoc RootBase::shortname
    std::string shortnameImpl() const
    {
      return classname();
    }

    /// \copydoc RootBase::createData
    ConstraintData createDataImpl() const
    {
      return ConstraintData(*this);
    }

    // Sizes -------------------------

    /// \copydoc RootBase::maxResidualSizeImpl
    template<ConstraintSelectionType Sel>
    int residualSizeImpl(ConstraintSelectionTag<Sel> sel) const
    {
      PINOCCHIO_UNUSED_VARIABLE(sel);
      return int(m_active_dofs.size());
    }

    /// \copydoc RootBase::symmetricConeResidualSize
    template<ConstraintSelectionType Sel>
    int symmetricConeResidualSizeImpl(ConstraintSelectionTag<Sel> sel) const
    {
      PINOCCHIO_UNUSED_VARIABLE(sel);
      return 2 * residualSize();
    }

    /// \copydoc RootBase::symmetricConeResidualScalingSize
    template<ConstraintSelectionType Sel>
    int symmetricConeResidualScalingSizeImpl(ConstraintSelectionTag<Sel> sel) const
    {
      PINOCCHIO_UNUSED_VARIABLE(sel);
      return 2 * residualSize();
    }

    // Methods for algorithms -------------

    /// \copydoc RootBase::set
    ConstraintSet setImpl(const ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return ConstraintSet(m_friction_lower_limit, m_friction_upper_limit);
    }

    /// \copydoc RootBase::calc
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
    }

    /// \copydoc RootBase::getRowSparsityPattern
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void getRowSparsityPatternImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id,
      BooleanVector & result) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < residualSize());
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);

      result = model.sparsity_pattern_vector[m_active_joint_ids[size_t(row_id)]];
    }

    /// \copydoc RootBase::getRowIndexes
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void getRowIndexesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id,
      EigenIndexVector & result) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < residualSize());
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);

      result = model.span_indexes_vector[m_active_joint_ids[size_t(row_id)]];
    }

    /// \copydoc RootBase::jacobian
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename JacobianMatrix>
    void jacobianImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const;

    /// \copydoc RootBase::jacobianMatrixProduct
    template<
      int OtherOptions,
      typename InputMatrix,
      template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(residualSize(), mat.cols());
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    /// \copydoc RootBase::jacobianMatrixProduct
    template<
      int OtherOptions,
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op>
    void jacobianMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot) const;

    /// \copydoc RootBase::jacobianTransposeMatrixProduct
    template<
      int OtherOptions,
      typename InputMatrix,
      template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef
        typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
          ReturnType;
      ReturnType res(model.nv, mat.cols());
      jacobianTransposeMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    /// \copydoc RootBase::jacobianTransposeMatrixProduct
    template<
      int OtherOptions,
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op>
    void jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot) const;

    /// \copydoc Base::mapConstraintForcesToJointTorques
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForcesLike,
      typename JointTorquesLike>
    void mapConstraintForceToJointTorquesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques) const;

    /// \copydoc Base::mapJointMotionsToConstraintMotions
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename JointMotionsLike,
      typename ConstraintMotionsLike>
    void mapJointMotionsToConstraintMotionImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions) const;

    /// \copydoc RootBase::appendCouplingConstraintInertias
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertiasImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const;

    /// \copydoc RootBase::appendCouplingConstraintInertias
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename MatrixOrMap,
      typename MapEnable,
      ReferenceFrame rf>
    void appendCouplingConstraintInertiasImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const MatrixBlockElementTpl<MatrixOrMap, MapEnable> & constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const;

  protected:
    // ------------------------------
    // PROTECTED METHODS
    // ------------------------------

    /// \brief Initialization of the model.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void init(const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model);

    // ------------------------------
    // MEMBERS
    // ------------------------------

    JointIndexVector m_active_joints;
    EigenIndexVector m_active_dofs;
    JointIndexVector m_active_joint_ids;

    VectorXs m_friction_lower_limit;
    VectorXs m_friction_upper_limit;

    using BaseCommonParameters::m_compliance;
  }; // struct JointFrictionConstraintModelTpl

  template<typename _Scalar, int _Options>
  struct JointFrictionConstraintDataTpl
  : ConstraintDataBase<JointFrictionConstraintDataTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Type defs
    // --------------------------------------------------------------
    // CRTP related types -------------------------------------------
    typedef JointFrictionConstraintDataTpl Self;
    typedef ConstraintDataBase<Self> Base;

    // Retrieving traits --------------------------------------------
    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintData ConstraintData;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    // Base usage ---------------------------------------------------
    using Base::classname;

    // Useful types ------------------------------------------------
    typedef std::vector<JointIndex> JointIndexVector;

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

    // Constructors ------------------

    /// \brief Default constructor
    JointFrictionConstraintDataTpl()
    {
    }

    /// \brief Constructor from a constraint model
    explicit JointFrictionConstraintDataTpl(const ConstraintModel & /*cmodel*/)
    {
    }

    // Operators ---------------------

    /// \brief Comparison operator
    bool operator==(const JointFrictionConstraintDataTpl & /*other*/) const
    {
      return true;
    }

    /// \brief Comparison operator
    bool operator!=(const JointFrictionConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    // General -----------------------

    /// \copydoc Base::classname
    static std::string classnameImpl()
    {
      return std::string("JointFrictionConstraintData");
    }

    /// \copydoc Base::shortname
    std::string shortnameImpl() const
    {
      return classname();
    }
  }; // struct JointFrictionConstraintDataTpl

  template<typename Scalar, int Options>
  template<int OtherOptions, template<typename, int> class JointCollectionTpl>
  void JointFrictionConstraintModelTpl<Scalar, Options>::init(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
  {
    m_active_dofs.reserve(size_t(model.nv));
    for (const JointIndex joint_id : m_active_joints)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        joint_id < model.joints.size(),
        "joint_id is larger than the total number of joints contained in the model.");

      const auto nv = model.nvs[joint_id];
      const auto idx_v = model.idx_vs[joint_id];

      for (int k = 0; k < nv; ++k)
      {
        const int row_id = idx_v + k;
        m_active_dofs.push_back(row_id);
        m_active_joint_ids.push_back(joint_id);
      }
    }

    {
      // Fill lower/upper bound based on input model
      // These values can be changed later if needed
      m_friction_lower_limit.resize(static_cast<Eigen::Index>(m_active_dofs.size()));
      m_friction_upper_limit.resize(static_cast<Eigen::Index>(m_active_dofs.size()));
      Eigen::Index idx = 0;
      for (const auto dof : m_active_dofs)
      {
        m_friction_lower_limit.coeffRef(idx) = model.lowerDryFrictionLimit.coeff(dof);
        m_friction_upper_limit.coeffRef(idx) = model.upperDryFrictionLimit.coeff(dof);
        ++idx;
      }
    }

    m_compliance = ResidualVectorType::Zero(residualSize());
  }

  // -------------------------------
  // IMPLEMENTATIONS OF BASE METHODS
  // -------------------------------

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    template<typename, int> class JointCollectionTpl,
    typename JacobianMatrix>
  void JointFrictionConstraintModelTpl<Scalar, Options>::jacobianImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & /*data*/,
    const ConstraintData & /*cdata*/,
    const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
  {
    JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

    const JointFrictionConstraintModelTpl & cmodel = *this;
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.rows(), cmodel.residualSize(),
      "The input/output Jacobian matrix does not have the right number of rows.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      jacobian_matrix.cols(), model.nv,
      "The input/output Jacobian matrix does not have the right number of cols.");

    jacobian_matrix.setZero();
    for (size_t row_id = 0; row_id < m_active_dofs.size(); ++row_id)
    {
      const auto col_id = m_active_dofs[row_id];
      jacobian_matrix(Eigen::Index(row_id), col_id) = Scalar(1);
    }
  }

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointFrictionConstraintModelTpl<Scalar, Options>::jacobianMatrixProductImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), residualSize());
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    for (size_t row_id = 0; row_id < m_active_dofs.size(); ++row_id)
    {
      const auto col_id = m_active_dofs[row_id];

      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(Eigen::Index(row_id)) -= mat.row(col_id);
      else
        res.row(Eigen::Index(row_id)) += mat.row(col_id);
    }
  }

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    typename InputMatrix,
    typename OutputMatrix,
    template<typename, int> class JointCollectionTpl,
    AssignmentOperatorType op>
  void JointFrictionConstraintModelTpl<Scalar, Options>::jacobianTransposeMatrixProductImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<InputMatrix> & mat,
    const Eigen::MatrixBase<OutputMatrix> & _res,
    AssignmentOperatorTag<op> aot) const
  {
    OutputMatrix & res = _res.const_cast_derived();

    PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), residualSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(aot);

    if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
      res.setZero();

    for (size_t row_id = 0; row_id < m_active_dofs.size(); ++row_id)
    {
      const auto col_id = m_active_dofs[row_id];

      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
        res.row(col_id) -= mat.row(Eigen::Index(row_id));
      else
        res.row(col_id) += mat.row(Eigen::Index(row_id));
    }
  }

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintForcesLike,
    typename JointTorquesLike>
  void JointFrictionConstraintModelTpl<Scalar, Options>::mapConstraintForceToJointTorquesImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
    const Eigen::MatrixBase<JointTorquesLike> & joint_torques_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), residualSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_torques_.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);

    auto & joint_torques = joint_torques_.const_cast_derived();

    for (size_t dof_id = 0; dof_id < m_active_dofs.size(); ++dof_id)
    {
      const auto row_id = m_active_dofs[dof_id];

      joint_torques.row(row_id) += constraint_forces.row(Eigen::Index(dof_id));
    }
  }

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    template<typename, int> class JointCollectionTpl,
    typename JointMotionsLike,
    typename ConstraintMotionsLike>
  void JointFrictionConstraintModelTpl<Scalar, Options>::mapJointMotionsToConstraintMotionImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
    const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions_) const
  {
    PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions_.rows(), residualSize());
    PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.rows(), model.nv);
    PINOCCHIO_UNUSED_VARIABLE(data);
    PINOCCHIO_UNUSED_VARIABLE(cdata);

    auto & constraint_motions = constraint_motions_.const_cast_derived();

    for (size_t dof_id = 0; dof_id < m_active_dofs.size(); ++dof_id)
    {
      const auto row_id = m_active_dofs[dof_id];

      constraint_motions.row(Eigen::Index(dof_id)) = joint_motions.row(row_id);
    }
  }

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    template<typename, int> class JointCollectionTpl,
    typename VectorNLike,
    ReferenceFrame rf>
  void JointFrictionConstraintModelTpl<Scalar, Options>::appendCouplingConstraintInertiasImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
    const ReferenceFrameTag<rf> reference_frame) const
  {
    PINOCCHIO_UNUSED_VARIABLE(cdata);
    PINOCCHIO_UNUSED_VARIABLE(reference_frame);

    PINOCCHIO_CHECK_ARGUMENT_SIZE(
      diagonal_constraint_inertia.size(), residualSize(),
      "The diagonal_constraint_inertia is of wrong size.");

    Eigen::Index row_id = 0;
    for (const JointIndex joint_id : m_active_joints)
    {
      const auto joint_nv = model.nvs[joint_id];
      const auto joint_diagonal_constraint_inertia =
        diagonal_constraint_inertia.segment(row_id, joint_nv);

      data.joint_apparent_inertia[joint_id].diagonal() += joint_diagonal_constraint_inertia;

      row_id += joint_nv;
    }
  }

  template<typename Scalar, int Options>
  template<
    int OtherOptions,
    template<typename, int> class JointCollectionTpl,
    typename MatrixOrMap,
    typename MapEnable,
    ReferenceFrame rf>
  void JointFrictionConstraintModelTpl<Scalar, Options>::appendCouplingConstraintInertiasImpl(
    const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
    DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
    const ConstraintData & cdata,
    const MatrixBlockElementTpl<MatrixOrMap, MapEnable> & constraint_inertia,
    const ReferenceFrameTag<rf> reference_frame) const
  {
    assert(constraint_inertia.size() == residualSize());
    switch (constraint_inertia.type())
    {
    case MatrixBlockType::Zero: {
      break;
    }
    case MatrixBlockType::Identity: {
      appendCouplingConstraintInertiasImpl(
        model, data, cdata, VectorXs::Ones(residualSize()), reference_frame);
      break;
    }
    case MatrixBlockType::ScalarIdentity: {
      const Scalar val = constraint_inertia.container()(0, 0);
      appendCouplingConstraintInertiasImpl(
        model, data, cdata, VectorXs::Constant(residualSize(), val), reference_frame);
      break;
    }
    case MatrixBlockType::Diagonal: {
      appendCouplingConstraintInertiasImpl(
        model, data, cdata, constraint_inertia.container().col(0), reference_frame);
      break;
    }
    case MatrixBlockType::Plain: {
      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument,
        "JointFrictionConstraintModel does not support Plain inertia blocks.");
      break;
    }
    default:
      assert(false && "Invalid MatrixBlockType for JointFrictionConstraintModel.");
      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument, "Invalid MatrixBlockType for JointFrictionConstraintModel.");
    }
  }

} // namespace pinocchio
