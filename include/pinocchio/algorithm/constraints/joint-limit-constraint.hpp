//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
#define __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__

#include "pinocchio/math/fwd.hpp"

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/orthant-cone.hpp"

#include "pinocchio/algorithm/constraints/jointwise-constraint-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

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
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::POSITION_LEVEL;
    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_baumgarte_corrector_vector = false;
    static constexpr bool constant_size = false;

    typedef JointLimitConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef NonNegativeOrthantConeTpl<Scalar> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef VectorXs VectorConstraintSize;

    typedef VectorXs ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef EigenStorageTpl<VectorXs> EigenStorageVector;
    typedef typename EigenStorageVector::RefMapType ActiveComplianceVectorTypeRef;
    typedef typename EigenStorageVector::ConstRefMapType ActiveComplianceVectorTypeConstRef;

    typedef VectorXs BaumgarteVectorType;
    typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorType>
      BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersRef;
    typedef const BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersConstRef;

    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::
        Matrix<Scalar, Eigen::Dynamic, InputMatrix::ColsAtCompileTime, InputMatrixPlain::Options>
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
    typedef typename traits<Self>::EigenStorageVector EigenStorageVector;
    typedef typename traits<Self>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef
      typename traits<Self>::ActiveComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;
    typedef
      typename traits<Self>::BaumgarteCorrectorVectorParameters BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::ConstraintSet ConstraintSet;

    using RootBase::jacobian;
    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    typedef std::vector<BooleanVector> VectorOfBooleanVector;
    typedef std::vector<EigenIndexVector> VectofOfEigenIndexVector;
    typedef std::vector<size_t> VectorOfSize;
    typedef std::vector<JointIndex> JointIndexVector;
    typedef typename traits<Self>::VectorXs VectorXs;
    typedef typename traits<Self>::RowVectorXs RowVectorXs;
    typedef VectorXs VectorConstraintSize;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      CompactTangentMap;

    JointLimitConstraintModelTpl()
    {
    }

    JointLimitConstraintModelTpl(const JointLimitConstraintModelTpl & other)
    {
      *this = other;
    }

    template<template<typename, int> class JointCollectionTpl>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & _activable_joints)
    : joint_nqs(model.nqs)
    , joint_nvs(model.nvs)
    , joint_idx_vs(model.idx_vs)
    {
      init(
        model, _activable_joints, model.lowerPositionLimit, model.upperPositionLimit,
        model.positionLimitMargin);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & _activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub)
    : joint_nqs(model.nqs)
    , joint_nvs(model.nvs)
    , joint_idx_vs(model.idx_vs)
    {
      init(model, _activable_joints, lb, ub, model.positionLimitMargin);
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration,
      typename VectorMarginConfiguration>
    JointLimitConstraintModelTpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & _activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
      const Eigen::MatrixBase<VectorMarginConfiguration> & margin)
    : joint_nqs(model.nqs)
    , joint_nvs(model.nvs)
    , joint_idx_vs(model.idx_vs)
    {
      init(model, _activable_joints, lb, ub, margin);
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, JointLimitConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, JointLimitConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.activable_joints = activable_joints;
      res.nq_reduce = nq_reduce;
      res.nv_max_atom = nv_max_atom;
      res.lower_activable_size = lower_activable_size;
      res.row_sparsity_pattern = row_sparsity_pattern;
      res.row_indexes = row_indexes;
      res.position_limit = position_limit.template cast<NewScalar>();
      res.position_margin = position_margin.template cast<NewScalar>();
      res.activable_idx_qs = activable_idx_qs;
      res.activable_idx_rows = activable_idx_rows;
      res.activable_idx_qs_reduce = activable_idx_qs_reduce;
      res.joint_nqs = joint_nqs;
      res.joint_nvs = joint_nvs;
      res.joint_idx_vs = joint_idx_vs;

      return res;
    }

    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    int activeSize(const ConstraintData & constraint_data) const
    {
      return int(constraint_data.active_idx_rows.size());
    }

    int size() const
    {
      return int(activable_idx_rows.size());
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    BaseCommonParameters & base_common_parameters()
    {
      return static_cast<BaseCommonParameters &>(*this);
    }
    const BaseCommonParameters & base_common_parameters() const
    {
      return static_cast<const BaseCommonParameters &>(*this);
    }

    /// \brief Compute the constraint residual for the active constraints.
    /// \note `calc` internally calls `resize`.
    /// \note the constraint residual is computed based on the model's lower/upper position limits,
    /// joint limit margin. It is also computed based on data.q_in.
    template<template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

    /// \brief Resize the constraint if needed at the current state given by data and store the
    /// results in cdata.
    template<template<typename, int> class JointCollectionTpl>
    void resize(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const;

    /// \brief Returns the sparsity associated with a given row
    const BooleanVector & getRowSparsityPattern(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_sparsity_pattern[activable_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the sparsity associated with a given row
    /// This vector is computed when calling the calc method.
    const BooleanVector & getActiveRowSparsityPattern(
      const ConstraintData & constraint_data, const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < activeSize(constraint_data));
      return row_sparsity_pattern[constraint_data.active_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the vector of the activable indexes associated with a given row
    const EigenIndexVector & getActivableRowIndexes(const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < size());
      return row_indexes[activable_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    /// This vector is computed when calling the calc method.
    const EigenIndexVector & getActiveRowIndexes(
      const ConstraintData & constraint_data, const Eigen::DenseIndex row_id) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id < activeSize(constraint_data));
      return row_indexes[constraint_data.active_idx_rows[static_cast<size_t>(row_id)]];
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeConstRef
    getActivecompliance(const ConstraintData & constraint_data) const
    {
      return constraint_data.active_compliance;
    }

    /// \brief Returns the compliance internally stored in the constraint model
    ActiveComplianceVectorTypeRef getActivecompliance(ConstraintData & constraint_data) const
    {
      return constraint_data.active_compliance;
    }

    /// \copydoc Base::set
    ConstraintSet setImpl() const
    {
      return ConstraintSet();
    }

    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const;

    /// \copydoc Base::mapConstraintForcesToJointTorques
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForcesLike,
      typename JointTorquesLike>
    void mapConstraintForceToJointTorques(
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
    void mapJointMotionsToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions) const;

    /// \brief Returns the vector of the active indexes associated with a given row
    /// This vector is computed when calling the calc method.
    const VectorOfSize & getActiveSetIndexes(const ConstraintData & constraint_data) const
    {
      return constraint_data.active_set_indexes;
    }

    /// Specialized accessor

    const JointIndexVector & getActivableJoints() const
    {
      return activable_joints;
    }
    int getNqReduce() const
    {
      return nq_reduce;
    }
    int getNvMaxAtom() const
    {
      return nv_max_atom;
    }
    int lowerSize() const
    {
      return lower_activable_size;
    }
    int lowerActiveSize(const ConstraintData & constraint_data) const
    {
      return constraint_data.lower_active_size;
    }
    int upperSize() const
    {
      return size() - lowerSize();
    }
    int upperActiveSize(const ConstraintData & constraint_data) const
    {
      return activeSize(constraint_data) - lowerActiveSize(constraint_data);
    }

    const VectorXs & getBoundPositionLimit() const
    {
      return position_limit;
    }
    const VectorXs & getBoundPositionMargin() const
    {
      return position_margin;
    }

    const EigenIndexVector & getActivableIdxQs() const
    {
      return activable_idx_qs;
    }

    const EigenIndexVector & getActivableIdxQsReduce() const
    {
      return activable_idx_qs_reduce;
    }

    const EigenIndexVector & getActiveIdxQsReduce(const ConstraintData & constraint_data) const
    {
      return constraint_data.active_idx_qs_reduce;
    }

    // row_sparsity_pattern, row_indexes, activable_idx_rows, active_idx_rows are
    // not exposed as they only privately allow getRowActiv[e/able]SparsityPattern and
    // getRowActiv[e/able]Indexes

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other JointLimitConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const JointLimitConstraintModelTpl & other) const
    {
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && activable_joints == other.activable_joints && nq_reduce == other.nq_reduce
             && nv_max_atom == other.nv_max_atom
             && lower_activable_size == other.lower_activable_size
             && row_sparsity_pattern == other.row_sparsity_pattern
             && row_indexes == other.row_indexes && position_limit == other.position_limit
             && position_margin == other.position_margin
             && activable_idx_qs == other.activable_idx_qs
             && activable_idx_rows == other.activable_idx_rows
             && activable_idx_qs_reduce == other.activable_idx_qs_reduce
             && joint_nqs == other.joint_nqs && joint_nvs == other.joint_nvs
             && joint_idx_vs == other.joint_idx_vs;
    }

    bool operator!=(const JointLimitConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    JointLimitConstraintModelTpl & operator=(const JointLimitConstraintModelTpl & other)
    {
      if (this != &other)
      {
        base_common_parameters() = other.base_common_parameters();
        activable_joints = other.activable_joints;
        nq_reduce = other.nq_reduce;
        nv_max_atom = other.nv_max_atom;
        lower_activable_size = other.lower_activable_size;
        row_sparsity_pattern = other.row_sparsity_pattern;
        row_indexes = other.row_indexes;
        position_limit = other.position_limit;
        position_margin = other.position_margin;
        activable_idx_qs = other.activable_idx_qs;
        activable_idx_rows = other.activable_idx_rows;
        activable_idx_qs_reduce = other.activable_idx_qs_reduce;
        joint_nqs = other.joint_nqs;
        joint_nvs = other.joint_nvs;
        joint_idx_vs = other.joint_idx_vs;
      }
      return *this;
    }

    // Jacobian operations

    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const;

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(size(), mat.cols());
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const;

    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
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

    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const;

  public:
    static std::string classname()
    {
      return std::string("JointLimitConstraintModel");
    }
    std::string shortname() const
    {
      return classname();
    }

  protected:
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorLowerConfiguration,
      typename VectorUpperConfiguration,
      typename VectorMarginConfiguration>
    void init(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const JointIndexVector & _activable_joints,
      const Eigen::MatrixBase<VectorLowerConfiguration> & lb,
      const Eigen::MatrixBase<VectorUpperConfiguration> & ub,
      const Eigen::MatrixBase<VectorMarginConfiguration> & margin);

    /// @brief List of activable joints, i.e. joint can that can be indeed reach its bounds. size
    /// nja
    JointIndexVector activable_joints;

    /// @brief nq size given the considered joints
    /// nq_reduce = SUM(j in activable_joints) j.nq
    int nq_reduce;

    /// @brief maximum nv size of all atomic (even if there is some composite) in all
    /// activable_joints
    int nv_max_atom;

    /// @brief number of activable lower bound limits
    int lower_activable_size;

    /// @brief Sparsity pattern for each considered joint. size nja
    VectorOfBooleanVector row_sparsity_pattern;
    VectofOfEigenIndexVector row_indexes;

    /// @brief Limit value of lower and upper bound in the constraint (size size()=lsize+usize)
    VectorXs position_limit;
    /// @brief Margin value of lower and upper bound in the constraint (size size()=lsize+usize)
    VectorXs position_margin;

    /// @brief give for each activable constraint the index in the configuration vector
    EigenIndexVector activable_idx_qs;

    /// @brief give for each activable constraint the row_id of sparsity pattern
    VectorOfSize activable_idx_rows;

    /// @brief give for each activable constraint of sparsity pattern
    EigenIndexVector activable_idx_qs_reduce;

    std::vector<int> joint_nqs, joint_nvs, joint_idx_vs;

    using BaseCommonParameters::m_baumgarte_parameters;
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

    JointLimitConstraintDataTpl()
    : constraint_residual(constraint_residual_storage.map())
    , active_compliance(active_compliance_storage.map())
    {
    }

    JointLimitConstraintDataTpl(const JointLimitConstraintDataTpl & other)
    : constraint_residual(constraint_residual_storage.map())
    , active_compliance(active_compliance_storage.map())
    {
      *this = other;
    }

    explicit JointLimitConstraintDataTpl(const ConstraintModel & constraint_model)
    : compact_tangent_map(
        CompactTangentMap::Zero(constraint_model.getNqReduce(), constraint_model.getNvMaxAtom()))
    , rowise_tangent_map(
        size_t(constraint_model.getNqReduce()), size_t(constraint_model.getNvMaxAtom()))
    , activable_constraint_residual(constraint_model.size())
    , constraint_residual_storage(constraint_model.size())
    , constraint_residual(constraint_residual_storage.map())
    , active_compliance(active_compliance_storage.map())
    {

      // Allocate the maximum size for the dynamic quantities
      lower_active_size = 0;
      const auto max_size = size_t(constraint_model.activable_idx_rows.size());

      active_set_indexes.reserve(max_size);
      active_idx_rows.reserve(max_size);
      active_idx_qs_reduce.reserve(max_size);
      active_compliance_storage.reserve(int(max_size));
      assert(
        constraint_model.activeSize(*this) == constraint_model.lowerActiveSize(*this)
        == constraint_model.upperActiveSize(*this) == 0);

      constraint_residual_storage.resize(0);

      // Allocate slices for rowise_tangent_map
      for (const auto joint_id : constraint_model.activable_joints)
      {
        const auto joint_nq = constraint_model.joint_nqs[joint_id];
        const auto joint_nv = constraint_model.joint_nvs[joint_id];
        for (int i = 0; i < joint_nq; ++i)
          rowise_tangent_map.push_back(1, joint_nv);
      }

      assert(rowise_tangent_map.size() == size_t(constraint_model.getNqReduce()));
    }

    bool operator==(const JointLimitConstraintDataTpl & other) const
    {
      if (this == &other)
        return true;
      return (
        rowise_tangent_map == other.rowise_tangent_map
        && activable_constraint_residual == other.activable_constraint_residual
        && constraint_residual_storage == other.constraint_residual_storage
        && constraint_residual == other.constraint_residual

        && lower_active_size == other.lower_active_size
        && active_set_indexes == other.active_set_indexes
        && active_idx_rows == other.active_idx_rows
        && active_idx_qs_reduce == other.active_idx_qs_reduce

        && active_compliance_storage == other.active_compliance_storage
        && active_compliance == other.active_compliance);
    }

    bool operator!=(const JointLimitConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    JointLimitConstraintDataTpl & operator=(const JointLimitConstraintDataTpl & other)
    {
      if (this != &other)
      {
        compact_tangent_map = other.compact_tangent_map;
        rowise_tangent_map = other.rowise_tangent_map;
        activable_constraint_residual = other.activable_constraint_residual;
        constraint_residual_storage = other.constraint_residual_storage;
        constraint_residual = constraint_residual_storage.map();

        lower_active_size = other.lower_active_size;
        active_set_indexes = other.active_set_indexes;
        active_idx_rows = other.active_idx_rows;
        active_idx_qs_reduce = other.active_idx_qs_reduce;

        active_compliance_storage = other.active_compliance_storage;
        active_compliance = active_compliance_storage.map();
      }
      return *this;
    }

    /// @brief Compact storage of the tangent map
    CompactTangentMap compact_tangent_map;

    RowVectorStack rowise_tangent_map;

    /// \brief Residual of all the activable constraints
    VectorXs activable_constraint_residual;

    /// \brief Residual of the active constraints
    EigenStorageVector constraint_residual_storage;
    typename EigenStorageVector::RefMapType constraint_residual;

  protected:
    /// @brief number of active lower bound limits activable
    int lower_active_size;

    /// \brief Vector containing the indexes of the constraints in the active set.
    /// the size of the vector is cmodel.size(cdata)
    /// each element have value < size()
    VectorOfSize active_set_indexes;

    /// @brief give for each active constraint the row_id of sparsity pattern
    VectorOfSize active_idx_rows;

    /// @brief give for each active constraint of sparsity pattern
    EigenIndexVector active_idx_qs_reduce;

    /// \brief Compliance of the active constraints
    EigenStorageVector active_compliance_storage;
    typename EigenStorageVector::RefMapType active_compliance;

  public:
    static std::string classname()
    {
      return std::string("JointLimitConstraintData");
    }
    std::string shortname() const
    {
      return classname();
    }
  };
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/joint-limit-constraint.hxx"

#endif // ifndef __pinocchio_algorithm_constraints_joint_limit_constraint_hpp__
