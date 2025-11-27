//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
#define __pinocchio_algorithm_constraint_model_generic_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"
#include "pinocchio/algorithm/constraints/visitors/constraint-model-visitor.hpp"

namespace pinocchio
{

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct traits<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
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

    // static constexpr ConstraintFormulationLevel NOT USED;
    static constexpr ConstraintBehaviour constraint_behaviour = ConstraintBehaviour::GENERAL;

    // The generic behave as if it has those elements but raise an error if underlying type do not
    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_compliance_member = true;
    static constexpr bool has_set = true;

    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintModel;
    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef boost::blank ConstraintSet;
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
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef VectorXs ComplianceVectorType;
    typedef Eigen::Ref<ComplianceVectorType> ComplianceVectorTypeRef;
    typedef Eigen::Ref<const ComplianceVectorType> ComplianceVectorTypeConstRef;

    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    // Will be removed
    typedef ComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef ComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;

    // Not needed anymore
    // typedef VectorXs BaumgarteVectorType;
    // typedef Eigen::Ref<VectorXs> BaumgarteVectorTypeRef;
    // typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorTypeRef>
    //   BaumgarteCorrectorVectorParameters;
    // typedef BaumgarteCorrectorVectorParameters BaumgarteCorrectorVectorParametersRef;
    // typedef Eigen::Ref<const VectorXs> BaumgarteVectorTypeConstRef;
    // typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorTypeConstRef>
    //   BaumgarteCorrectorVectorParametersConstRef;
  };

  template<
    typename _Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl>
  struct ConstraintModelTpl
  : ConstraintModelBase<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  , ConstraintCollectionTpl<_Scalar, _Options>::ConstraintModelVariant
  , serialization::Serializable<ConstraintModelTpl<_Scalar, _Options, ConstraintCollectionTpl>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef ConstraintModelTpl Self;
    typedef ConstraintModelBase<Self> Base;
    typedef ConstraintModelBase<Self> RootBase;

    typedef ConstraintDataTpl<Scalar, Options, ConstraintCollectionTpl> ConstraintData;
    typedef ConstraintCollectionTpl<Scalar, Options> ConstraintCollection;
    typedef typename ConstraintCollection::ConstraintDataVariant ConstraintDataVariant;
    typedef typename ConstraintCollection::ConstraintModelVariant ConstraintModelVariant;
    typedef typename traits<Self>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Self>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef
      typename traits<Self>::ActiveComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;
    typedef typename Base::BaumgarteCorrectorParameters BaumgarteCorrectorParameters;

    using typename Base::BooleanVector;
    using typename Base::EigenIndexVector;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    ConstraintModelTpl()
    : ConstraintModelVariant()
    {
    }

    ConstraintModelTpl(const ConstraintModelVariant & cmodel_variant)
    : ConstraintModelVariant(cmodel_variant)
    {
    }

    template<typename ContraintModelDerived>
    ConstraintModelTpl(const ConstraintModelBase<ContraintModelDerived> & cmodel)
    : ConstraintModelVariant((ConstraintModelVariant)cmodel.derived())
    {
      BOOST_MPL_ASSERT(
        (boost::mpl::contains<typename ConstraintModelVariant::types, ContraintModelDerived>));
    }

    ConstraintModelVariant & toVariant()
    {
      return static_cast<ConstraintModelVariant &>(*this);
    }

    const ConstraintModelVariant & toVariant() const
    {
      return static_cast<const ConstraintModelVariant &>(*this);
    }

    template<typename ConstraintModelDerived>
    bool isEqual(const ConstraintModelBase<ConstraintModelDerived> & other) const
    {
      return ::pinocchio::isEqual(*this, other.derived());
    }

    bool isEqual(const ConstraintModelTpl & other) const
    {
      return toVariant() == other.toVariant();
    }

    /// \brief Comparison operator
    bool operator==(const ConstraintModelTpl & other) const
    {
      return isEqual(other);
    }

    /// \brief Comparison operator
    bool operator!=(const ConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc RootBase::classname
    static std::string classnameImpl()
    {
      return "ConstraintModel";
    }

    /// \copydoc RootBase::shortname
    std::string shortnameImpl() const
    {
      return ::pinocchio::visitors::shortname(*this);
    }

    /// \copydoc RootBase::size
    int sizeImpl() const
    {
      return ::pinocchio::visitors::size(*this);
    }

    /// \copydoc RootBase::activeSize
    int activeSizeImpl(const ConstraintData & constraint_data) const
    {
      return ::pinocchio::visitors::activeSize(*this, constraint_data);
    }

    /// \copydoc RootBase::set
    boost::blank setImpl() const
    {
      static boost::blank val;
      PINOCCHIO_THROW_PRETTY(
        std::runtime_error, "Set method is not accessible for ConstraintModelTpl.");
      return val;
    }

    /// \copydoc RootBase::createData
    ConstraintData createDataImpl() const
    {
      return ::pinocchio::visitors::createData<Scalar, Options, ConstraintCollectionTpl>(*this);
    }

    /// \copydoc RootBase::calc
    template<template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      ::pinocchio::visitors::calc(*this, model, data, cdata);
    }

    /// \copydoc RootBase::jacobian
    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobianImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix) const
    {
      ::pinocchio::visitors::jacobian(
        *this, model, data, cdata, jacobian_matrix.const_cast_derived());
    }

    /// \copydoc RootBase::getActivableRowIndexes
    const EigenIndexVector & getActivableRowIndexesImpl(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getActivableRowIndexes(*this, row_id);
    }

    /// \copydoc RootBase::getActiveRowIndexes
    const EigenIndexVector & getActiveRowIndexesImpl(
      const ConstraintData & constraint_data, const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getActiveRowIndexes(*this, constraint_data, row_id);
    }

    /// \copydoc RootBase::getRowSparsityPattern
    const BooleanVector & getRowSparsityPatternImpl(const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getRowSparsityPattern(*this, row_id);
    }

    /// \copydoc RootBase::getActiveRowSparsityPattern
    const BooleanVector & getActiveRowSparsityPatternImpl(
      const ConstraintData & constraint_data, const Eigen::DenseIndex row_id) const
    {
      return ::pinocchio::visitors::getActiveRowSparsityPattern(*this, constraint_data, row_id);
    }

    /// \copydoc RootBase::getActiveCompliance
    ActiveComplianceVectorTypeConstRef getActivecomplianceImpl(const ConstraintData & cdata) const
    {
      return ::pinocchio::visitors::getActiveCompliance(*this, cdata);
    }

    /// \copydoc RootBase::getActiveCompliance
    ActiveComplianceVectorTypeRef getActivecomplianceImpl(ConstraintData & cdata) const
    {
      return ::pinocchio::visitors::getActiveCompliance(*this, cdata);
    }

    /// \copydoc RootBase::jacobianMatrixProduct
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(activeSize(cdata), input_matrix.cols());
      jacobianMatrixProduct(model, data, cdata, input_matrix.derived(), res);
      return res;
    }

    /// \copydoc RootBase::jacobianMatrixProduct
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      ::pinocchio::visitors::jacobianMatrixProduct(
        *this, model, data, cdata, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
    }

    /// \copydoc RootBase::jacobianTransposeMatrixProduct
    template<
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & input_matrix) const
    {
      typedef
        typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
          ReturnType;
      ReturnType res(model.nv, input_matrix.cols());
      jacobianTransposeMatrixProduct(*this, model, data, cdata, input_matrix.derived(), res);
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
      const Eigen::MatrixBase<InputMatrix> & input_matrix,
      const Eigen::MatrixBase<OutputMatrix> & result_matrix,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      ::pinocchio::visitors::jacobianTransposeMatrixProduct(
        *this, model, data, cdata, input_matrix.derived(), result_matrix.const_cast_derived(), aot);
    }

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
      const ReferenceFrameTag<rf> reference_frame) const
    {
      ::pinocchio::visitors::appendCouplingConstraintInertias(
        *this, model, data, cdata, diagonal_constraint_inertia.derived(), reference_frame);
    }

    /// \copydoc RootBase::mapConstraintForceToJointSpace
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForceLike,
      typename ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    void mapConstraintForceToJointSpaceImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
      ReferenceFrameTag<rf> reference_frame) const
    {
      ::pinocchio::visitors::mapConstraintForceToJointSpace(
        *this, model, data, cdata, constraint_forces, joint_forces,
        joint_torques.const_cast_derived(), reference_frame);
    }

    /// \copydoc RootBase::mapJointSpaceToConstraintMotion
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename JointMotionsLike,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointSpaceToConstraintMotionImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<JointMotionsLike> & joint_generalized_velocity,
      const Eigen::MatrixBase<VectorLike> & constraint_motions,
      ReferenceFrameTag<rf> reference_frame) const
    {
      ::pinocchio::visitors::mapJointSpaceToConstraintMotion(
        *this, model, data, cdata, joint_motions, joint_generalized_velocity,
        constraint_motions.const_cast_derived(), reference_frame);
    }

    /// \copydoc RootBase::compliance_impl
    ComplianceVectorTypeConstRef compliance_impl() const
    {
      return ::pinocchio::visitors::compliance(*this);
    }

    /// \copydoc RootBase::compliance_impl
    ComplianceVectorTypeRef compliance_impl()
    {
      return ::pinocchio::visitors::compliance(*this);
    }

    // CHOICE: right now we use the scalar Baumgarte
    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersConstRef baumgarte_corrector_vector_parameters_impl() const
    // {
    //   return ::pinocchio::visitors::getBaumgarteCorrectorVectorParameters(*this);
    // }

    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersRef baumgarte_corrector_vector_parameters_impl()
    // {
    //   return ::pinocchio::visitors::getBaumgarteCorrectorVectorParameters(*this);
    // }

    /// \copydoc RootBase::baumgarte_corrector_parameters
    const BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl() const
    {
      return ::pinocchio::visitors::getBaumgarteCorrectorParameters(*this);
    }

    /// \copydoc RootBase::baumgarte_corrector_parameters
    BaumgarteCorrectorParameters & baumgarte_corrector_parameters_impl()
    {
      return ::pinocchio::visitors::getBaumgarteCorrectorParameters(*this);
    }

  }; // struct ConstraintModelTpl
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_model_generic_hpp__
