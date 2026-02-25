//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/jointwise-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"
#include "pinocchio/algorithm/constraints/sets/orthant-cone.hpp"
#include "pinocchio/algorithm/constraints/sets/orthant-cone-jordan-operation.hpp"

#include "pinocchio/container/eigen-storage.hpp"
#include "pinocchio/container/matrix-stack.hpp"

namespace pinocchio
{

  // --------------------------------------------------------------
  // Cast
  // --------------------------------------------------------------
  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, JointLimitConstraintModelTpl<Scalar, Options>>
  {
    typedef JointLimitConstraintModelTpl<NewScalar, Options> type;
  };

  // --------------------------------------------------------------
  // Traits
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef JointLimitConstraintModelTpl<_Scalar, _Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<_Scalar, _Options> ConstraintData;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    // --------------------------------------------------------------
    // Traits characterizing the constraints
    // --------------------------------------------------------------
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::POSITION_LEVEL;
    static constexpr ConstraintSizeType constraint_size_type = ConstraintSizeType::BOUNDED;

    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_set = true;
    static constexpr bool is_inequality_constraint = true;

    // --------------------------------------------------------------
    // Traits for associated struct and sizes
    // --------------------------------------------------------------
    typedef NonNegativeOrthantConeTpl<Scalar> ConstraintSet;
    typedef NonNegativeOrthantJordanOperationTpl<Scalar, Options> JordanOperation;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

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
  struct traits<JointLimitConstraintDataTpl<_Scalar, _Options>>
  : traits<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
  };

  // --------------------------------------------------------------
  // Unsafe wrapper
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct Unsafe<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    typedef JointLimitConstraintModelTpl<_Scalar, _Options> SafeSelf;
    typedef typename SafeSelf::VectorOfSize VectorOfSize;

    explicit Unsafe(SafeSelf & self)
    : self(self)
    {
    }

    // Non-const getter of active_idx_in_activable for custom selection
    VectorOfSize & active_idx_in_activable()
    {
      return self.m_cursel_active_idx_in_activable;
    }

  protected:
    SafeSelf & self;
  };

  // --------------------------------------------------------------
  // Struct
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct JointLimitConstraintModelTpl
  : JointWiseConstraintModelBase<JointLimitConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<JointLimitConstraintModelTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Type defs
    // --------------------------------------------------------------
    // CRTP related types -------------------------------------------
    typedef JointLimitConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelCommonParameters<Self> BaseCommonParameters;
    typedef ConstraintModelBase<Self> RootBase;

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

    typedef typename traits<Self>::RowVectorXs RowVectorXs;

    // Friendship ---------------------------------------------------
    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintModelTpl;

    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintDataTpl;

    friend struct Unsafe<Self>;

    // Base usage ---------------------------------------------------
    using RootBase::classname;
    using RootBase::jacobianMatrixProduct;
    using RootBase::jacobianTransposeMatrixProduct;
    using RootBase::residualSize;
    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    // Useful types ------------------------------------------------
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef EigenStorageTpl<VectorXs> EigenStorageVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      CompactTangentMap;
    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectorOfEigenIndexVector;
    typedef std::vector<size_t> VectorOfSize;
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

    // Unsafe API --------------------
    Unsafe<Self> unsafe()
    {
      return Unsafe<Self>(*this);
    }

    // Constructors ------------------

    /// \brief Default constructor
    JointLimitConstraintModelTpl()
    {
    }

    /// \brief Constructor from model only.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    {
      size_t n_joints = model.joints.size();
      JointIndexVector activable_joints;
      activable_joints.reserve(n_joints);
      for (size_t i = 1; i < n_joints; ++i)
      {
        activable_joints.push_back(static_cast<JointIndex>(i));
      }
      init(
        model, activable_joints, model.lowerPositionLimit, model.upperPositionLimit,
        model.positionLimitMargin);
    }

    /// \brief Constructor from model, activable joints, lower, upper and margin joint limits.
    /// \note lb, ub and margin must be of size nq. They are the bounds of the entire model.
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration,
      typename VectorMarginConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
      const Eigen::MatrixBase<VectorMarginConfiguration> & margin)
    : Base(model)
    {
      init(model, activable_joints, lb, ub, margin);
    }

    /// \brief Constructor from model, activable joints, lower and upper joint limits.
    /// \note lb and ub must be of size nq. They are the bounds of the entire model.
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
    : JointLimitConstraintModelTpl(model, activable_joints, lb, ub, model.positionLimitMargin)
    {
    }

    /// \brief Constructor from model and activable joints.
    /// Activable joints are joints that can become active/non-active
    /// depending on their position w.r.t the joint limit margin.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints)
    : JointLimitConstraintModelTpl(
        model,
        activable_joints,
        model.lowerPositionLimit,
        model.upperPositionLimit,
        model.positionLimitMargin)
    {
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
      // Constraint definition related
      res.m_selected_joints = m_selected_joints;
      res.m_selected_row_sparsity_pattern = m_selected_row_sparsity_pattern;
      res.m_selected_row_indexes = m_selected_row_indexes;
      res.m_selected_joint_nqs = m_selected_joint_nqs;
      res.m_selected_joint_nvs = m_selected_joint_nvs;
      res.m_selected_joint_idx_vs = m_selected_joint_idx_vs;
      res.m_nq_reduce = m_nq_reduce;
      res.m_max_of_nvs = m_max_of_nvs;
      res.m_lower_activable_residual_size = m_lower_activable_residual_size;
      res.m_activable_idx_in_selected = m_activable_idx_in_selected;
      res.m_activable_idx_qs = m_activable_idx_qs;
      res.m_activable_idx_qs_reduce = m_activable_idx_qs_reduce;
      res.m_activable_position_limit = m_activable_position_limit.template cast<NewScalar>();
      res.m_activable_position_margin = m_activable_position_margin.template cast<NewScalar>();
      // Selection related
      res.m_cursel_active_idx_in_activable = m_cursel_active_idx_in_activable;
      res.m_cursel_lower_active_residual_size = m_cursel_lower_active_residual_size;
      res.m_cursel_active_idx_in_selected = m_cursel_active_idx_in_selected;
      res.m_cursel_active_idx_qs = m_cursel_active_idx_qs;
      res.m_cursel_active_idx_qs_reduce = m_cursel_active_idx_qs_reduce;

      return res;
    }

    ///
    /// \brief Comparison operator
    ///
    /// \param[in] other Other JointLimitConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributes
    /// must be the same).
    ///
    bool operator==(const JointLimitConstraintModelTpl & other) const
    {
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && m_selected_joints == other.m_selected_joints
             && m_selected_row_sparsity_pattern == other.m_selected_row_sparsity_pattern
             && m_selected_row_indexes == other.m_selected_row_indexes
             && m_selected_joint_nqs == other.m_selected_joint_nqs
             && m_selected_joint_nvs == other.m_selected_joint_nvs
             && m_selected_joint_idx_vs == other.m_selected_joint_idx_vs
             && m_nq_reduce == other.m_nq_reduce && m_max_of_nvs == other.m_max_of_nvs
             && m_lower_activable_residual_size == other.m_lower_activable_residual_size
             && m_activable_idx_in_selected == other.m_activable_idx_in_selected
             && m_activable_idx_qs == other.m_activable_idx_qs
             && m_activable_idx_qs_reduce == other.m_activable_idx_qs_reduce
             && m_activable_position_limit == other.m_activable_position_limit
             && m_activable_position_margin == other.m_activable_position_margin
             && m_cursel_active_idx_in_activable == other.m_cursel_active_idx_in_activable
             && m_cursel_lower_active_residual_size == other.m_cursel_lower_active_residual_size
             && m_cursel_active_idx_in_selected == other.m_cursel_active_idx_in_selected
             && m_cursel_active_idx_qs == other.m_cursel_active_idx_qs
             && m_cursel_active_idx_qs_reduce == other.m_cursel_active_idx_qs_reduce;
    }

    /// \brief Comparison operator
    bool operator!=(const JointLimitConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    /// Specialized accessors for constraints definition --------

    /// \copydoc m_selected_joints
    const JointIndexVector & getSelectedJoints() const
    {
      return m_selected_joints;
    }
    /// \copydoc m_nq_reduce
    int getNqReduce() const
    {
      return m_nq_reduce;
    }
    /// \copydoc m_max_of_nvs
    int getMaxOfNvs() const
    {
      return m_max_of_nvs;
    }
    /// \copydoc m_activable_position_limit
    const VectorXs & getActivablePositionLimit() const
    {
      return m_activable_position_limit;
    }
    /// \copydoc m_activable_position_margin
    const VectorXs & getActivablePositionMargin() const
    {
      return m_activable_position_margin;
    }

    // m_selected_row_sparsity_pattern, m_selected_row_indexes,
    // m_selected_joint_nqs, m_selected_joint_nvs, m_selected_joint_idx_vs,
    // m_activable_idx_in_selected, m_activable_idx_qs, m_activable_idx_qs_reduce
    // m_cursel_active_idx_in_selected, m_cursel_active_idx_qs, m_cursel_active_idx_qs_reduce
    // not exposed as they are only used privately.

    /// Specialized methods for lower and upper sizes

    /// \brief Return the maximum residual size of constraints that are lower limits
    template<ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
    int lowerResidualSize(ConstraintSelectionTag<Sel> sel = CurrentSelection()) const
    {
      PINOCCHIO_UNUSED_VARIABLE(sel);
      if constexpr (std::is_same_v<ConstraintSelectionTag<Sel>, MaximalSelection>)
      {
        return m_lower_activable_residual_size;
      }
      else // Current selection
      {
        return m_cursel_lower_active_residual_size;
      }
    }

    /// \brief Return the maximum residual size of constraints that are upper limits
    template<ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
    int upperResidualSize(ConstraintSelectionTag<Sel> sel = CurrentSelection()) const
    {
      return residualSize(sel) - lowerResidualSize(sel);
    }

    /// Specialized methods for constraints definition --------

    /// \brief Set activable_[position_limit|margin] of size residualSize(MaximalSelection) from lb,
    /// ub, margin of size model.nq
    /// \note Expect a limit or margin vector of size model.nq
    template<
      typename VectorLike1,
      typename VectorLike2,
      typename VectorLike3,
      ConstraintSelectionType Sel = ConstraintSelectionType::MAXIMAL>
    void setPositionLimitAndMargin(
      const Eigen::MatrixBase<VectorLike1> & lb,
      const Eigen::MatrixBase<VectorLike2> & ub,
      const Eigen::MatrixBase<VectorLike3> & margin,
      ConstraintSelectionTag<Sel> sel = MaximalSelection());

    /// \brief Reset the current selection to the maximal selection.
    ///
    /// This method sets the active constraint selection to include all activable constraints,
    /// effectively making the current selection equal to the maximal selection. After calling
    /// this method, residualSize(CurrentSelection()) == residualSize(MaximalSelection()).
    ///
    /// \note This is useful when you want to consider all possible joint limit constraints
    /// without any filtering based on the current configuration.
    void makeSelectionMaximal();

    ///
    /// \brief Update the current selection to include only constraints that are near their limits.
    ///
    /// This method selects constraints where the joint configuration is within the margin distance
    /// from either the lower or upper position limit. Constraints outside this margin are excluded
    /// from the current selection.
    ///
    /// For lower bounds: activates if q[idx] - lower_limit <= margin
    /// For upper bounds: activates if upper_limit - q[idx] <= margin
    ///
    /// \param[in] q The joint configuration vector of size model.nq used to evaluate proximity
    ///              to joint limits.
    ///
    /// \note This method updates the current selection state (m_cursel_* members) and affects
    ///       the result of residualSize(CurrentSelection()).
    /// \note After calling this method, only constraints near their limits will be active.
    /// \note To restore all constraints, call makeSelectionMaximal().
    ///
    template<typename VectorLike>
    void makeSelectionFilteredByLimitProximity(const Eigen::MatrixBase<VectorLike> & q);

    const VectorOfSize & active_idx_in_activable() const
    {
      return m_cursel_active_idx_in_activable;
    }

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

    // Sizes -------------------------

    /// \copydoc RootBase::residualSizeImpl
    template<ConstraintSelectionType Sel>
    int residualSizeImpl(ConstraintSelectionTag<Sel> sel) const
    {
      PINOCCHIO_UNUSED_VARIABLE(sel);
      if constexpr (std::is_same_v<ConstraintSelectionTag<Sel>, MaximalSelection>)
      {
        return int(m_activable_idx_in_selected.size());
      }
      else // CurrentSelection
      {
        return int(m_cursel_active_idx_in_selected.size());
      }
    }

    /// \copydoc RootBase::symmetricConeResidualSize
    template<ConstraintSelectionType Sel>
    int symmetricConeResidualSizeImpl(ConstraintSelectionTag<Sel> sel) const
    {
      return residualSize(sel);
    }

    /// \copydoc RootBase::symmetricConeResidualScalingSize
    template<ConstraintSelectionType Sel>
    int symmetricConeResidualScalingSizeImpl(ConstraintSelectionTag<Sel> sel) const
    {
      return residualSize(sel);
    }

    // Hyperparameters handling -----------

    /// \brief Set the compliance
    template<typename VectorLike, ConstraintSelectionType Sel>
    void
    setComplianceImpl(const Eigen::MatrixBase<VectorLike> & vector, ConstraintSelectionTag<Sel> sel)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(vector.size()) == residualSize(sel));
      if constexpr (std::is_same_v<ConstraintSelectionTag<Sel>, MaximalSelection>)
      {
        m_compliance = vector;
      }
      else // Current selection
      {
        for (Eigen::Index row_id = 0; row_id < residualSize(); ++row_id)
        {
          const Eigen::Index idx =
            Eigen::Index(m_cursel_active_idx_in_activable[static_cast<size_t>(row_id)]);
          m_compliance[idx] = vector[row_id];
        }
      }
    }

    /// \copydoc RootBase::retrieveCompliance
    template<typename VectorLike, ConstraintSelectionType Sel>
    void retrieveComplianceImpl(
      const Eigen::MatrixBase<VectorLike> & res_, ConstraintSelectionTag<Sel> sel) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(res_.size()) == residualSize(sel));
      auto res = res_.const_cast_derived();
      if constexpr (std::is_same_v<ConstraintSelectionTag<Sel>, MaximalSelection>)
      {
        res = m_compliance;
      }
      else // Current selection
      {
        for (Eigen::Index row_id = 0; row_id < residualSize(); ++row_id)
        {
          const Eigen::Index idx =
            Eigen::Index(m_cursel_active_idx_in_activable[static_cast<size_t>(row_id)]);
          res[row_id] = m_compliance[idx];
        }
      }
    }

    // Methods for algorithms -------------

    /// \copydoc RootBase::set
    ConstraintSet setImpl(const ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return ConstraintSet();
    }

    /// \copydoc RootBase::calc
    /// \note the constraint residual is computed based on the model's lower/upper position limits,
    /// joint limit margin and data.q_in.
    /// \note it calls computeResidualAndSelectActiveConstraints which select the constraints to
    /// consider
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

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

    /// \copydoc RootBase::getRowSparsityPattern
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    const BooleanVector & getRowSparsityPatternImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(row_id) < residualSize());
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      const size_t idx = m_cursel_active_idx_in_selected[static_cast<size_t>(row_id)];
      return m_selected_row_sparsity_pattern[idx];
    }

    /// \copydoc RootBase::getRowIndexes
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    const EigenIndexVector & getRowIndexesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(int(row_id) < residualSize());
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      const size_t idx = m_cursel_active_idx_in_selected[static_cast<size_t>(row_id)];
      return m_selected_row_indexes[idx];
    }

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
      const std::vector<MatrixBlockElementTpl<MatrixOrMap, MapEnable>> & constraint_inertias,
      const ReferenceFrameTag<rf> reference_frame,
      std::size_t & inner_constraint_id) const;

  protected:
    // ------------------------------
    // PROTECTED METHODS
    // ------------------------------

    /// \brief Initialize the constraint model with model, activable joints, lower, upper and margin
    /// of joint limits.
    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration,
      typename VectorMarginConfiguration>
    void init(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndexVector & activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
      const Eigen::MatrixBase<VectorMarginConfiguration> & margin);

    // ------------------------------
    // MEMBERS
    // ------------------------------

    // Maximal selection definition ---------------------------------------------------

    /// \brief List of selected joints, i.e. joint can that can be indeed reach its bounds.
    /// size = nActivableJoints
    JointIndexVector m_selected_joints;

    /// \brief Sparsity pattern for each selected joints.
    /// size = nActivableJoints
    VectorOfBooleanVector m_selected_row_sparsity_pattern;
    VectorOfEigenIndexVector m_selected_row_indexes;

    /// \brief Vector of size for selected joints
    std::vector<int> m_selected_joint_nqs, m_selected_joint_nvs, m_selected_joint_idx_vs;

    /// \brief nq size given the selected joints
    /// m_nq_reduce = SUM(j in m_selected_joints) j.nq
    int m_nq_reduce;

    /// \brief max of nv of the selected joints
    /// m_max_of_nvs = MAX(j in m_selected_joints) j.nv
    int m_max_of_nvs;

    /// \brief number of activable lower bound limits. By convention for i=0..lmrs lower limit and
    /// i=lmrs..mrs upper limits.
    // m_lower_activable_residual_size <= residualSize(MaximalSelection())
    int m_lower_activable_residual_size;

    /// \brief give for each activable constraint the index of related joint in m_selected_joints
    /// size = residualSize(MaximalSelection())
    VectorOfSize m_activable_idx_in_selected;

    /// \brief give for each activable constraint the index in [0, Nq] for the activable constraint
    /// size = residualSize(MaximalSelection())
    EigenIndexVector m_activable_idx_qs;

    /// \brief give for each activable constraint the index in [0, Nqred] for the activable
    /// size = residualSize(MaximalSelection())
    EigenIndexVector m_activable_idx_qs_reduce;

    /// \brief Limit value of lower and upper bound in the constraint (size size()=lsize+usize)
    /// size = residualSize(MaximalSelection())
    VectorXs m_activable_position_limit;

    /// \brief Margin value of lower and upper bound in the constraint (size size()=lsize+usize)
    /// size = residualSize(MaximalSelection())
    VectorXs m_activable_position_margin;

    /// \brief Baumgarte correction parameters of the constraint model
    using BaseCommonParameters::m_baumgarte_parameters;

    /// \brief Compliance of the constraint model
    /// size = residualSize(MaximalSelection())
    using BaseCommonParameters::m_compliance;

    // Current selection definition ---------------------------------------------------

    /// \brief Vector containing the indexes of the activable constraints that are currently used.
    /// The size of the vector denoted, residualSize, is given by
    /// cmodel.residualSize(CurrentSelection()) each element have value <
    /// cmodel.residualSize(MaximalSelection()) This vector totally define the state of the
    /// constraint size = residualSize(CurrentSelection())
    VectorOfSize m_cursel_active_idx_in_activable;

    /// \brief number of active lower bound limits activable (<= residualSize)
    /// It is the number of element in active_idx_in_activable that are <
    /// m_lower_activable_residual_size
    // m_cursel_lower_active_residual_size <= m_lower_activable_residual_size
    // m_cursel_lower_active_residual_size <= residualSize(CurrentSelection())
    int m_cursel_lower_active_residual_size;

    /// \brief Proxys to avoid calculus in algorithmic methods
    /// size = residualSize(CurrentSelection())
    VectorOfSize m_cursel_active_idx_in_selected;
    EigenIndexVector m_cursel_active_idx_qs;
    EigenIndexVector m_cursel_active_idx_qs_reduce;
  }; // struct JointLimitConstraintModelTpl

  template<typename _Scalar, int _Options>
  struct JointLimitConstraintDataTpl
  : ConstraintDataBase<JointLimitConstraintDataTpl<_Scalar, _Options>>
  {
    // --------------------------------------------------------------
    // Type defs
    // --------------------------------------------------------------
    // CRTP related types -------------------------------------------
    typedef JointLimitConstraintDataTpl Self;
    typedef ConstraintDataBase<Self> Base;

    // Retrieving traits --------------------------------------------
    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintData ConstraintData;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    typedef typename ConstraintModel::VectorXs VectorXs;
    typedef typename ConstraintModel::RowVectorXs RowVectorXs;
    typedef typename ConstraintModel::CompactTangentMap CompactTangentMap;
    typedef typename ConstraintModel::EigenStorageVector EigenStorageVector;
    typedef typename ConstraintModel::BooleanVector BooleanVector;
    typedef typename ConstraintModel::EigenIndexVector EigenIndexVector;
    typedef typename ConstraintModel::VectorOfSize VectorOfSize;
    typedef typename ConstraintModel::JointIndexVector JointIndexVector;

    // Friendship ---------------------------------------------------
    template<typename NewScalar, int NewOptions>
    friend struct JointLimitConstraintModelTpl;

    // Base usage ---------------------------------------------------
    using Base::classname;

    // Useful types ------------------------------------------------
    typedef MatrixStackTpl<RowVectorXs> RowVectorStack;

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
    {
    }

    /// \brief Copy constructor
    JointLimitConstraintDataTpl(const JointLimitConstraintDataTpl & other)
    : constraint_residual_storage(other.constraint_residual_storage)
    , compact_tangent_map(other.compact_tangent_map)
    , rowise_tangent_map(other.rowise_tangent_map)
    {
    }

    /// \brief Move constructor
    JointLimitConstraintDataTpl(JointLimitConstraintDataTpl && other)
    : constraint_residual_storage(other.constraint_residual_storage)
    , compact_tangent_map(std::move(other.compact_tangent_map))
    , rowise_tangent_map(other.rowise_tangent_map)
    {
    }

    /// \brief Constructor from a constraint model
    explicit JointLimitConstraintDataTpl(const ConstraintModel & cmodel)
    : constraint_residual_storage(cmodel.residualSize(MaximalSelection()))
    , compact_tangent_map(CompactTangentMap::Zero(cmodel.m_nq_reduce, cmodel.m_max_of_nvs))
    , rowise_tangent_map(
        static_cast<size_t>(cmodel.m_nq_reduce), static_cast<size_t>(cmodel.m_max_of_nvs))
    {
      // Allocate slices for rowise_tangent_map
      for (size_t sel_id = 0; sel_id < cmodel.m_selected_joints.size(); sel_id++)
      {
        const int joint_nq = cmodel.m_selected_joint_nqs[sel_id];
        const int joint_nv = cmodel.m_selected_joint_nvs[sel_id];
        for (int i = 0; i < joint_nq; ++i)
          rowise_tangent_map.push_back(1, joint_nv);
      }

      assert(rowise_tangent_map.size() == static_cast<size_t>(cmodel.m_nq_reduce));
    }

    // Operators ---------------------

    /// Custom copy operator because of the ref. Thus we apply rule of five
    /// \brief Copy operator
    JointLimitConstraintDataTpl & operator=(const JointLimitConstraintDataTpl & other)
    {
      if (this != &other)
      {
        constraint_residual_storage = other.constraint_residual_storage;
        compact_tangent_map = other.compact_tangent_map;
        rowise_tangent_map = other.rowise_tangent_map;
      }
      return *this;
    }

    /// \brief Move assignment operator
    JointLimitConstraintDataTpl & operator=(JointLimitConstraintDataTpl && other)
    {
      if (this != &other)
      {
        constraint_residual_storage =
          other.constraint_residual_storage; // Move not defined for EigenStorage
        compact_tangent_map = std::move(other.compact_tangent_map);
        rowise_tangent_map = other.rowise_tangent_map; // Move not defined for MatrixStack
      }
      return *this;
    }

    /// \brief Default destructor
    ~JointLimitConstraintDataTpl() = default;

    /// \brief Comparison operator
    bool operator==(const JointLimitConstraintDataTpl & other) const
    {
      if (this == &other)
        return true;
      return (
        constraint_residual_storage == other.constraint_residual_storage
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

    /// \brief Residual of the active constraints
    /// size = cmodel.residualSize(CurrentSelection())
    /// capacity = cmodel.residualSize(MaximalSelection())
    EigenStorageVector constraint_residual_storage;
    typename EigenStorageVector::RefMapType constraint_residual = constraint_residual_storage.map();

    /// \brief Compact storages of the tangent map
    CompactTangentMap compact_tangent_map;
    RowVectorStack rowise_tangent_map;
  }; // struct JointLimitConstraintDataTpl

} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-limit-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
