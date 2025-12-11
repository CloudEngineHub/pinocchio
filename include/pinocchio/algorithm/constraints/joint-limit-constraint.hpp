//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/jointwise-constraint-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"
#include "pinocchio/algorithm/constraints/sets/orthant-cone.hpp"

#include "pinocchio/container/eigen-storage.hpp"
#include "pinocchio/container/matrix-stack.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, JointLimitConstraintModelTpl<Scalar, Options>>
  {
    typedef JointLimitConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
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
      ConstraintFormulationLevel::POSITION_LEVEL;
    static constexpr ConstraintSizeType constraint_size_type = ConstraintSizeType::BOUNDED;

    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_compliance_member = true;
    static constexpr bool has_set = true;

    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef NonNegativeOrthantConeTpl<Scalar> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    // --------------------------------------------------------------
    // Traits for the algorithmic methods on current state
    // --------------------------------------------------------------
    // Elementary types
    typedef Eigen::Matrix<Scalar, Size, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Eigen::Matrix<Scalar, Size, 1, Options> VectorConstraintSize;

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

    // -------------------------------
    // Traits for holded Data
    // -------------------------------
    typedef Eigen::Matrix<Scalar, Size, 1, Options> ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointLimitConstraintDataTpl<_Scalar, _Options>>
  : traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
  };

  template<typename _Scalar, int _Options>
  struct JointLimitConstraintModelTpl
  : JointWiseConstraintModelBase<JointLimitConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef JointLimitConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelBase<Self> RootBase;

    typedef ConstraintModelCommonParameters<JointLimitConstraintModelTpl> BaseCommonParameters;

    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintDataTpl;

    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintModelTpl;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<JointLimitConstraintModelTpl>::constraint_formulation_level;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::ConstraintSet ConstraintSet;

    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    typedef typename traits<Self>::VectorXs VectorXs;
    typedef typename traits<Self>::RowVectorXs RowVectorXs;

    typedef EigenStorageTpl<VectorXs> EigenStorageVector;
    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;
    typedef std::vector<size_t> VectorOfSize;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef VectorXs VectorConstraintSize;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      CompactTangentMap;

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
    JointLimitConstraintModelTpl()
    {
    }

    /// \brief Copy constructor
    JointLimitConstraintModelTpl(const JointLimitConstraintModelTpl & other)
    {
      *this = other;
    }

    /// \brief Constructor from model and activable joints.
    /// Activable joints are joints that can become active/non-active
    /// depending on their position w.r.t the joint limit margin.
    template<template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints)
    {
      init(
        model, activable_joints, model.lowerPositionLimit, model.upperPositionLimit,
        model.positionLimitMargin);
    }

    /// \brief Constructor from model, activable joints, lower and upper joint limits.
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
    {
      init(model, activable_joints, lb, ub, model.positionLimitMargin);
    }

    /// \brief Constructor from model, activable joints, lower, upper and margin joint limits.
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration,
      typename VectorMarginConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
      const Eigen::MatrixBase<VectorMarginConfiguration> & margin)
    {
      init(model, activable_joints, lb, ub, margin);
    }

    // Operators ---------------------

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, JointLimitConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, JointLimitConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.selected_joints = selected_joints;
      res.selected_row_sparsity_pattern = selected_row_sparsity_pattern;
      res.selected_row_indexes = selected_row_indexes;
      res.selected_joint_nqs = selected_joint_nqs;
      res.selected_joint_nvs = selected_joint_nvs;
      res.selected_joint_idx_vs = selected_joint_idx_vs;
      res.nq_reduce = nq_reduce;
      res.max_of_nvs = max_of_nvs;
      res.activable_idx_in_selected = activable_idx_in_selected;
      res.activable_idx_qs = activable_idx_qs;
      res.activable_idx_qs_reduce = activable_idx_qs_reduce;
      res.activable_position_limit = activable_position_limit.template cast<NewScalar>();
      res.activable_position_margin = activable_position_margin.template cast<NewScalar>();
      res.lower_max_residual_size = lower_max_residual_size;

      return res;
    }

    ///
    /// \brief Comparison operator
    ///
    /// \param[in] other Other JointLimitConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const JointLimitConstraintModelTpl & other) const
    {
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && selected_joints == other.selected_joints
             && selected_row_sparsity_pattern == other.selected_row_sparsity_pattern
             && selected_row_indexes == other.selected_row_indexes
             && selected_joint_nqs == other.selected_joint_nqs
             && selected_joint_nvs == other.selected_joint_nvs
             && selected_joint_idx_vs == other.selected_joint_idx_vs && nq_reduce == other.nq_reduce
             && max_of_nvs == other.max_of_nvs
             && activable_idx_in_selected == other.activable_idx_in_selected
             && activable_idx_qs == other.activable_idx_qs
             && activable_idx_qs_reduce == other.activable_idx_qs_reduce
             && activable_position_limit == other.activable_position_limit
             && activable_position_margin == other.activable_position_margin
             && lower_max_residual_size == other.lower_max_residual_size;
    }

    /// \brief Comparison operator
    bool operator!=(const JointLimitConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    /// \brief Copy operator
    JointLimitConstraintModelTpl & operator=(const JointLimitConstraintModelTpl & other)
    {
      if (this != &other)
      {
        base_common_parameters() = other.base_common_parameters();
        selected_joints = other.selected_joints;
        selected_row_sparsity_pattern = other.selected_row_sparsity_pattern;
        selected_row_indexes = other.selected_row_indexes;
        selected_joint_nqs = other.selected_joint_nqs;
        selected_joint_nvs = other.selected_joint_nvs;
        selected_joint_idx_vs = other.selected_joint_idx_vs;
        nq_reduce = other.nq_reduce;
        max_of_nvs = other.max_of_nvs;
        activable_idx_in_selected = other.activable_idx_in_selected;
        activable_idx_qs = other.activable_idx_qs;
        activable_idx_qs_reduce = other.activable_idx_qs_reduce;
        activable_position_limit = other.activable_position_limit;
        activable_position_margin = other.activable_position_margin;
        lower_max_residual_size = other.lower_max_residual_size;
      }
      return *this;
    }

    /// Specialized accessors --------

    /// \copydoc selected_joints
    const JointIndexVector & getSelectedJoints() const
    {
      return selected_joints;
    }
    /// \copydoc nq_reduce
    int getNqReduce() const
    {
      return nq_reduce;
    }
    /// \copydoc max_of_nvs
    int getMaxOfNvs() const
    {
      return max_of_nvs;
    }
    /// \copydoc activable_position_limit
    const VectorXs & getActivablePositionLimit() const
    {
      return activable_position_limit;
    }
    /// \copydoc activable_position_margin
    const VectorXs & getActivablePositionMargin() const
    {
      return activable_position_margin;
    }
    /// \brief Return the maximum residual size of constraints that are lower limits
    int lowerMaxResidualSize() const
    {
      return lower_max_residual_size;
    }
    /// \brief Return the maximum residual size of constraints that are upper limits
    int upperMaxResidualSize() const
    {
      return maxResidualSize() - lowerMaxResidualSize();
    }
    /// \brief Return the residual size of constraints that are lower limits given the state of the
    /// constraint given by cdata
    int lowerResidualSize(const ConstraintData & cdata) const
    {
      return cdata.lower_residual_size;
    }
    /// \brief Return the residual size of constraints that are upper limits given the state of the
    /// constraint given by cdata
    int upperResidualSize(const ConstraintData & cdata) const
    {
      return residualSize(cdata) - lowerResidualSize(cdata);
    }

    VectorOfSize getActiveIdxInActivable(const ConstraintData & cdata) const
    {
      return cdata.active_idx_in_activable;
    }

    /// \brief Set activable_[position_limit|margin] of size maxResidualSize from lb, ub, margin of
    /// size model.nq
    /// \note Expect a limit or margin vector of size model.nq
    template<typename VectorLike1, typename VectorLike2, typename VectorLike3>
    void setPositionLimitAndMargin(
      const Eigen::MatrixBase<VectorLike1> & lb,
      const Eigen::MatrixBase<VectorLike2> & ub,
      const Eigen::MatrixBase<VectorLike3> & margin)
    {
      // Fill bound limit and margin for lower and upper for activable constraints
      activable_position_limit = VectorXs::Zero(Eigen::Index(maxResidualSize()));
      activable_position_margin = VectorXs::Zero(Eigen::Index(maxResidualSize()));
      Eigen::Index constraint_id = 0;
      for (; constraint_id < lowerMaxResidualSize(); ++constraint_id)
      {
        const Eigen::Index idx_q = activable_idx_qs[static_cast<size_t>(constraint_id)];
        activable_position_limit[constraint_id] = lb[idx_q];
        activable_position_margin[constraint_id] = margin[idx_q];
        assert(margin[idx_q] >= 0);
      }
      for (; constraint_id < maxResidualSize(); ++constraint_id)
      {
        const Eigen::Index idx_q = activable_idx_qs[static_cast<size_t>(constraint_id)];
        activable_position_limit[constraint_id] = ub[idx_q];
        activable_position_margin[constraint_id] = margin[idx_q];
        assert(margin[idx_q] >= 0);
      }
    }

    // selected_row_sparsity_pattern, selected_row_indexes,
    // activable_idx_in_selected, activable_idx_qs_reduce
    // not exposed as they only privately allow getRowActiv[e/able]SparsityPattern and
    // getRowActiv[e/able]Indexes

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    // General -----------------------

    /// \copydoc RootBase::classname
    static std::string classnameImpl()
    {
      return std::string("JointLimitConstraintModel");
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
      return int(activable_idx_in_selected.size());
    }

    // Methods for algorithms --------

    /// \copydoc RootBase::residualSize
    int residualSizeImpl(const ConstraintData & cdata) const
    {
      return int(cdata.active_idx_in_selected.size());
    }

    /// \copydoc RootBase::getRowSparsityPattern
    const BooleanVector &
    getRowSparsityPatternImpl(const ConstraintData & cdata, const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(row_id) < residualSize(cdata));
      const size_t idx = cdata.active_idx_in_selected[static_cast<size_t>(row_id)];
      return selected_row_sparsity_pattern[idx];
    }

    /// \copydoc RootBase::getRowIndexes
    const EigenIndexVector &
    getRowIndexesImpl(const ConstraintData & cdata, const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(row_id) < residualSize(cdata));
      const size_t idx = cdata.active_idx_in_selected[static_cast<size_t>(row_id)];
      return selected_row_indexes[idx];
    }

    /// \copydoc RootBase::set
    ConstraintSet setImpl() const
    {
      return ConstraintSet();
    }

    /// \copydoc RootBase::retrieveCompliance
    template<typename VectorLike>
    void retrieveComplianceImpl(
      const ConstraintData & cdata, const Eigen::MatrixBase<VectorLike> & res_) const
    {
      auto res = res_.const_cast_derived();
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(res.size()) == residualSize(cdata));
      for (Eigen::Index row_id = 0; row_id < residualSize(cdata); ++row_id)
      {
        const Eigen::Index idx =
          Eigen::Index(cdata.active_idx_in_activable[static_cast<size_t>(row_id)]);
        res[row_id] = m_compliance[idx];
      }
    }

    /// \copydoc RootBase::calc
    /// \note the constraint residual is computed based on the model's lower/upper position limits,
    /// joint limit margin and data.q_in.
    /// \note it calls computeResidualAndSelectActiveConstraints which select the constraints to
    /// consider
    template<template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

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
      ReturnType res(residualSize(cdata), mat.cols());
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

    /// \brief Initialize the constraint model with model, activable joints, lower, upper and margin
    /// of joint limits.
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration,
      typename VectorMarginConfiguration>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
      const Eigen::MatrixBase<VectorMarginConfiguration> & margin);

    /// \brief Resize the constraint if needed at the current state given by data and store the
    /// results in cdata.
    template<template<typename, int> class JointCollectionTpl>
    void computeResidualAndActiveConstraints(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

    // ------------------------------
    // MEMBERS
    // ------------------------------

    /// \brief List of selected joints, i.e. joint can that can be indeed reach its bounds.
    /// size = nActivableJoints
    JointIndexVector selected_joints;

    /// \brief Sparsity pattern for each selected joints.
    /// size = nActivableJoints
    VectorOfBooleanVector selected_row_sparsity_pattern;
    VectofOfEigenIndexVector selected_row_indexes;

    /// \brief Vector of size for selected joints
    std::vector<int> selected_joint_nqs, selected_joint_nvs, selected_joint_idx_vs;

    /// \brief nq size given the selected joints
    /// nq_reduce = SUM(j in selected_joints) j.nq
    int nq_reduce;

    /// \brief max of nv of the selected joints
    /// max_of_nvs = MAX(j in selected_joints) j.nv
    int max_of_nvs;

    /// \brief give for each activable constraint the index of related joint in selected_joints
    /// size = maxResidualSize
    VectorOfSize activable_idx_in_selected;

    /// \brief give for each activable constraint the index in [0, Nq] for the activable constraint
    /// size = maxResidualsize
    EigenIndexVector activable_idx_qs;

    /// \brief give for each activable constraint the index in [0, Nqred] for the activable
    /// constraint size = maxResidualsize
    EigenIndexVector activable_idx_qs_reduce;

    /// \brief Limit value of lower and upper bound in the constraint (size size()=lsize+usize)
    /// size = maxResidualsize
    VectorXs activable_position_limit;

    /// \brief Margin value of lower and upper bound in the constraint (size size()=lsize+usize)
    /// size = maxResidualsize
    VectorXs activable_position_margin;

    /// \brief number of activable lower bound limits. By convention for i=0..lmrs lower limit and
    /// i=lmrs..mrs upper limits. lower_max_residual_size <= maxResidualsize
    int lower_max_residual_size;

    /// \brief Baumgarte correction parameters of the constraint model
    using BaseCommonParameters::m_baumgarte_parameters;

    /// \brief Compliance of the constraint model
    /// size = maxResidualsize
    using BaseCommonParameters::m_compliance;
  };

  template<typename _Scalar, int _Options>
  struct JointLimitConstraintDataTpl
  : ConstraintDataBase<JointLimitConstraintDataTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef ConstraintDataBase<JointLimitConstraintDataTpl> Base;
    typedef std::vector<JointIndex> JointIndexVector;

    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintModelTpl;

    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;

    typedef typename ConstraintModel::VectorXs VectorXs;
    typedef typename ConstraintModel::RowVectorXs RowVectorXs;
    typedef typename ConstraintModel::CompactTangentMap CompactTangentMap;
    typedef typename ConstraintModel::EigenStorageVector EigenStorageVector;
    typedef typename ConstraintModel::BooleanVector BooleanVector;
    typedef typename ConstraintModel::EigenIndexVector EigenIndexVector;
    typedef typename ConstraintModel::VectorOfSize VectorOfSize;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;

    typedef MatrixStackTpl<RowVectorXs> RowVectorStack;

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
    JointLimitConstraintDataTpl()
    : constraint_residual(constraint_residual_storage.map())
    {
    }

    /// \brief Copy constructor
    JointLimitConstraintDataTpl(const JointLimitConstraintDataTpl & other)
    : constraint_residual(constraint_residual_storage.map())
    {
      *this = other;
    }

    /// \brief Constructor from a constraint_model
    explicit JointLimitConstraintDataTpl(const ConstraintModel & constraint_model)
    : activable_constraint_residual(constraint_model.maxResidualSize())
    , constraint_residual_storage(constraint_model.maxResidualSize())
    , constraint_residual(constraint_residual_storage.map())
    , compact_tangent_map(
        CompactTangentMap::Zero(constraint_model.getNqReduce(), constraint_model.getMaxOfNvs()))
    , rowise_tangent_map(
        static_cast<size_t>(constraint_model.getNqReduce()),
        static_cast<size_t>(constraint_model.getMaxOfNvs()))
    {
      // Allocate the maximum size for the dynamic quantities
      const size_t max_residual_size = static_cast<size_t>(constraint_model.maxResidualSize());
      active_idx_in_activable.reserve(max_residual_size);
      active_idx_in_selected.reserve(max_residual_size);
      active_idx_qs_reduce.reserve(max_residual_size);

      // By default all activable constraint are not active
      lower_residual_size = 0;
      constraint_residual_storage.resize(0);
      assert(
        constraint_model.residualSize(*this) == constraint_model.lowerResidualSize(*this)
        == constraint_model.upperResidualSize(*this) == 0);

      // Allocate slices for rowise_tangent_map
      for (size_t sel_id = 0; sel_id < constraint_model.selected_joints.size(); sel_id++)
      {
        const int joint_nq = constraint_model.selected_joint_nqs[sel_id];
        const int joint_nv = constraint_model.selected_joint_nvs[sel_id];
        for (int i = 0; i < joint_nq; ++i)
          rowise_tangent_map.push_back(1, joint_nv);
      }

      assert(rowise_tangent_map.size() == static_cast<size_t>(constraint_model.getNqReduce()));
    }

    // Operators ---------------------

    /// \brief Copy operator
    JointLimitConstraintDataTpl & operator=(const JointLimitConstraintDataTpl & other)
    {
      if (this != &other)
      {
        active_idx_in_activable = other.active_idx_in_activable;
        active_idx_in_selected = other.active_idx_in_selected;
        active_idx_qs_reduce = other.active_idx_qs_reduce;
        lower_residual_size = other.lower_residual_size;
        activable_constraint_residual = other.activable_constraint_residual;
        constraint_residual_storage = other.constraint_residual_storage;
        constraint_residual = constraint_residual_storage.map();
        compact_tangent_map = other.compact_tangent_map;
        rowise_tangent_map = other.rowise_tangent_map;
      }
      return *this;
    }

    /// \brief Comparison operator
    bool operator==(const JointLimitConstraintDataTpl & other) const
    {
      if (this == &other)
        return true;
      return (
        active_idx_in_activable == other.active_idx_in_activable
        && active_idx_in_selected == other.active_idx_in_selected
        && active_idx_qs_reduce == other.active_idx_qs_reduce
        && lower_residual_size == other.lower_residual_size
        && activable_constraint_residual == other.activable_constraint_residual
        && constraint_residual_storage == other.constraint_residual_storage
        && constraint_residual == other.constraint_residual
        && compact_tangent_map == other.compact_tangent_map
        && rowise_tangent_map == other.rowise_tangent_map);
    }

    /// \brief Comparison operator
    bool operator!=(const JointLimitConstraintDataTpl & other) const
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
      return std::string("JointLimitConstraintData");
    }

    /// \copydoc Base::shortname
    std::string shortnameImpl() const
    {
      return classname();
    }

    // ------------------------------
    // MEMBERS
    // ------------------------------
    // note: data is always public - use at your own risk

    /// \brief Vector containing the indexes of the activable constraints that are currently used.
    /// The size of the vector denoted, residualSize, is given by cmodel.residualSize(cdata)
    /// each element have value < cmodel.maxResidualSize()
    /// This vector totally define the state of the constraint
    /// size = residualSize
    VectorOfSize active_idx_in_activable;

    /// \brief give for each active constraint the row_id of sparsity pattern (size = csize)
    /// size = residualSize
    /// \note it is a proxy to avoid double derefrence as :
    /// activable_idx_in_selected[active_idx_in_activable[i]] = active_idx_in_selected[i]
    VectorOfSize active_idx_in_selected;

    /// \brief give for each active constraint of sparsity pattern (size = csize)
    /// size = residualSize
    /// \note it is a proxy to avoid double derefrence as :
    /// activable_idx_qs_reduce[active_idx_in_activable[i]] = active_idx_qs_reduce[i]
    EigenIndexVector active_idx_qs_reduce;

    /// \brief number of active lower bound limits activable (<= residualSize)
    /// It is the number of element in active_idx_in_activable that are < lower_max_residual_size
    int lower_residual_size;

    /// \brief Residual of all potential constraints
    /// size = maxResidualSize
    VectorXs activable_constraint_residual;

    /// \brief Residual of the active constraints
    /// size = residualSize
    /// capacity = maxResidualSize
    EigenStorageVector constraint_residual_storage;
    typename EigenStorageVector::RefMapType constraint_residual;

    /// \brief Compact storages of the tangent map
    CompactTangentMap compact_tangent_map;
    RowVectorStack rowise_tangent_map;
  };
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-limit-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
