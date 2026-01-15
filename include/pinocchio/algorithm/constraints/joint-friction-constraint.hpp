//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_friction_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_friction_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/jointwise-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/sets/box-set.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, JointFrictionConstraintModelTpl<Scalar, Options>>
  {
    typedef JointFrictionConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Traits characterizing the constraint behaviour in CRTP
    // --------------------------------------------------------------
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options,
      Size = Eigen::Dynamic
    };

    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::VELOCITY_LEVEL;
    static constexpr ConstraintSizeType constraint_size_type = ConstraintSizeType::CONSTANT;

    static constexpr bool has_baumgarte_corrector =
      false; // Baumgarte make sense and exist directly for the constraint
    static constexpr bool has_compliance_member =
      true; // The constraint itself posses a member m_compliance which can be set by the user
    static constexpr bool has_set = true; // The constraint itself defines the set, otherwise must
                                          // have a mechanism for set-related visitors

    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef JointFrictionConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointFrictionConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef BoxSetTpl<Scalar, Options> ConstraintSet;
    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    // --------------------------------------------------------------
    // Traits for the algorithmic methods on current state
    // --------------------------------------------------------------
    // Elementary types
    typedef Eigen::Matrix<Scalar, Size, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Eigen::Matrix<Scalar, Size, 1, Options> VectorConstraintSize;

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

    // -------------------------------
    // Traits for holded Data
    // -------------------------------
    typedef Eigen::Matrix<Scalar, Size, 1, Options> ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointFrictionConstraintDataTpl<_Scalar, _Options>>
  : traits<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  {
  };

  template<typename _Scalar, int _Options>
  struct JointFrictionConstraintModelTpl
  : JointWiseConstraintModelBase<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<JointFrictionConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    typedef JointFrictionConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelCommonParameters<Self> BaseCommonParameters;
    typedef ConstraintModelBase<JointFrictionConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct JointFrictionConstraintModelTpl;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<JointFrictionConstraintModelTpl>::constraint_formulation_level;
    typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::ConstraintSet ConstraintSet;

    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs VectorConstraintSize;

    using RootBase::classname;
    using RootBase::jacobianMatrixProduct;
    using RootBase::jacobianTransposeMatrixProduct;
    using RootBase::maxResidualSize;
    using RootBase::residualSize;

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

    /// \brief Constructor from model and m_active_joints.
    template<template<typename, int> class JointCollectionTpl>
    JointFrictionConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & m_active_joints)
    : m_active_joints(m_active_joints)
    {
      init(model, m_active_joints);
    }

    // Operators ---------------------

    ///
    /// \brief Comparison operator
    ///
    /// \param[in] other Other JointFrictionConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const JointFrictionConstraintModelTpl & other) const
    {
      if (this == &other)
        return true;
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && m_active_dofs == other.m_active_dofs
             && m_row_sparsity_pattern == other.m_row_sparsity_pattern
             && m_row_active_indexes == other.m_row_active_indexes
             && m_friction_lower_limit == other.m_friction_lower_limit
             && m_friction_upper_limit == other.m_friction_upper_limit;
    }

    /// \brief Comparison operator
    bool operator!=(const JointFrictionConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

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
      res.m_row_sparsity_pattern = m_row_sparsity_pattern;
      res.m_row_active_indexes = m_row_active_indexes;
      res.m_friction_lower_limit = m_friction_lower_limit.template cast<NewScalar>();
      res.m_friction_upper_limit = m_friction_upper_limit.template cast<NewScalar>();
      return res;
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
        lb.size() != maxResidualSize(), std::runtime_error, "lb should be the same as size()");
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
        ub.size() != maxResidualSize(), std::runtime_error, "ub should be the same as size()");
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

    // Size Management ---------------

    /// \copydoc RootBase::maxResidualSizeImpl
    int maxResidualSizeImpl() const
    {
      return int(m_active_dofs.size());
    }

    // Methods for algorithms --------

    /// \copydoc RootBase::set
    ConstraintSet setImpl(const ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return ConstraintSet(m_friction_lower_limit, m_friction_upper_limit);
    }

    /// \copydoc RootBase::getRowSparsityPattern
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

      return m_row_sparsity_pattern[size_t(row_id)];
    }

    /// \copydoc RootBase::getRowIndexes
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

      return m_row_active_indexes[size_t(row_id)];
    }

    /// \copydoc RootBase::calc
    template<template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
    }

    /// \copydoc RootBase::jacobian
    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobianImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const;

    /// \copydoc RootBase::jacobianMatrixProduct
    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(maxResidualSize(), mat.cols());
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    /// \copydoc RootBase::jacobianMatrixProduct
    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const;

    /// \copydoc RootBase::jacobianTransposeMatrixProduct
    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
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
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const;

    /// \copydoc Base::mapConstraintForcesToJointTorques
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForcesLike,
      typename JointTorquesLike>
    void mapConstraintForceToJointTorquesImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques) const;

    /// \copydoc Base::mapJointMotionsToConstraintMotions
    template<
      template<typename, int> class JointCollectionTpl,
      typename JointMotionsLike,
      typename ConstraintMotionsLike>
    void mapJointMotionsToConstraintMotionImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions) const;

    /// \copydoc RootBase::appendCouplingConstraintInertias
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertiasImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const;

  protected:
    // ------------------------------
    // PROTECTED METHODS
    // ------------------------------

    /// \brief Initialization of the model.
    template<template<typename, int> class JointCollectionTpl>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & m_active_joints);

    // ------------------------------
    // MEMBERS
    // ------------------------------

    JointIndexVector m_active_joints;
    EigenIndexVector m_active_dofs;
    VectorOfBooleanVector m_row_sparsity_pattern;
    VectofOfEigenIndexVector m_row_active_indexes;

    VectorXs m_friction_lower_limit;
    VectorXs m_friction_upper_limit;

    using BaseCommonParameters::m_compliance;
  };

  template<typename _Scalar, int _Options>
  struct JointFrictionConstraintDataTpl
  : ConstraintDataBase<JointFrictionConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    typedef ConstraintDataBase<JointFrictionConstraintDataTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;

    typedef JointFrictionConstraintModelTpl<Scalar, Options> ConstraintModel;

    using Base::classname;

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

    /// \brief Constructor from a constraint_model
    explicit JointFrictionConstraintDataTpl(const ConstraintModel & /*constraint_model*/)
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
  };

} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-friction-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_joint_friction_constraint_hpp__
