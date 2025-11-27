//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_base_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/common/model-entity.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"

namespace pinocchio
{

  template<typename Scalar>
  struct BaumgarteCorrectorParametersTpl;

  template<typename ConstraintDataDerived>
  struct ConstraintDataBase;

  enum struct ConstraintFormulationLevel
  {
    POSITION_LEVEL,    // scaling dt^2
    VELOCITY_LEVEL,    // scaling dt
    ACCELERATION_LEVEL // scaling 1
  };

  template<class Derived>
  struct ConstraintModelBase
  : NumericalBase<Derived>
  , ModelEntity<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options
    };
    static constexpr bool constant_size = traits<Derived>::constant_size;

    typedef ModelEntity<Derived> Base;

    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ConstraintSet ConstraintSet;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef typename traits<Derived>::ActiveComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef typename traits<Derived>::ActiveComplianceVectorTypeConstRef
      ActiveComplianceVectorTypeConstRef;

    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    // Not used anymore:
    // typedef typename traits<Derived>::BaumgarteCorrectorVectorParametersRef
    //   BaumgarteCorrectorVectorParametersRef;
    // typedef typename traits<Derived>::BaumgarteCorrectorVectorParametersConstRef
    //   BaumgarteCorrectorVectorParametersConstRef;

    // Usefull types ------------------------------------------------
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector;
    typedef std::vector<Eigen::DenseIndex> EigenIndexVector;

    using Base::createData;

    /// \brief Cast to derived class.
    Derived & derived()
    {
      return static_cast<Derived &>(*this);
    }

    /// \brief Const cast to derived class.
    const Derived & derived() const
    {
      return static_cast<const Derived &>(*this);
    }

    /// \brief Cast to NewScalar.
    template<typename NewScalar>
    typename CastType<NewScalar, Derived>::type cast() const
    {
      return derived().template cast<NewScalar>();
    }

    /// \brief Cast to NewScalar.
    template<typename OtherDerived>
    void cast(ConstraintModelBase<OtherDerived> & other) const
    {
      other.name = name;
    }

    /// \brief Cast to base.
    ConstraintModelBase & base()
    {
      return *this;
    }

    /// \brief Const cast to base.
    const ConstraintModelBase & base() const
    {
      return *this;
    }

    /// \brief Copy operator.
    template<typename OtherDerived>
    ConstraintModelBase & operator=(const ConstraintModelBase<OtherDerived> & other)
    {
      name = other.name;

      return *this;
    }

    /// \brief Equality comparison operator.
    template<typename OtherDerived>
    bool operator==(const ConstraintModelBase<OtherDerived> & other) const
    {
      return name == other.name;
    }

    /// \brief Difference comparison operator.
    template<typename OtherDerived>
    bool operator!=(const ConstraintModelBase<OtherDerived> & other) const
    {
      return !(*this == other);
    }

    /// \brief Returns a constraint data associated to this constraint model.
    ConstraintData createData() const
    {
      return derived().createDataImpl();
    }

    /// \brief Returns the name of the underlying class if this is a variant.
    std::string shortname() const
    {
      return derived().shortnameImpl();
    }

    /// \brief Returns the name of the underlying class if this is a variant.
    static std::string classname()
    {
      return Derived::classnameImpl();
    }

    /// \brief Prints the shortname of the constraint.
    void disp(std::ostream & os) const
    {
      using namespace std;
      os << shortname() << endl;
    }

    /// \copydoc disp
    friend std::ostream &
    operator<<(std::ostream & os, const ConstraintModelBase<Derived> & constraint)
    {
      constraint.disp(os);
      return os;
    }

    /// \brief Returns the (maximum) size of the constraint.
    int size() const
    {
      return derived().sizeImpl();
    }

    /// \brief Returns the current size of the constraint, typically after `calc` has been called.
    /// \note If constraints are dynamic (e.g. joint limits), activeSize is computed when
    /// calling the calc method.
    template<typename ConstraintDataDerived>
    int activeSize(const ConstraintDataBase<ConstraintDataDerived> & constraint_data) const
    {
      if constexpr (traits<Derived>::constant_size)
      {
        return size();
      }
      else
      {
        return derived().activeSizeImpl(constraint_data.derived());
      }
    }

    /// \brief Evaluate the constraint values at the current state given by data and store the
    /// results in cdata.
    template<int Options, template<typename, int> class JointCollectionTpl>
    void calc(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      derived().calcImpl(model, data, cdata);
    }

    /// \brief Evaluate the Jacobian associated to the constraint at the given state stored in data
    /// and cdata.
    /// The results Jacobian is evaluated in the jacobian input/output matrix.
    /// This method assumes that the constrained data is up-to-date (calc has been called).
    template<int Options, template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & jacobian_matrix) const
    {
      derived().jacobianImpl(model, data, cdata, jacobian_matrix.const_cast_derived());
    }

    /// \copydoc jacobian
    template<int Options, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::JacobianMatrixType jacobian(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      typedef typename traits<Derived>::JacobianMatrixType ReturnType;
      ReturnType res = ReturnType::Zero(activeSize(cdata), model.nv);

      jacobian(model, data, cdata, res);

      return res;
    }

    /// \brief Evaluates the constraint jacobian against an input matrix mat.
    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      return derived().jacobianMatrixProductImpl(model, data, cdata, mat.derived());
    }

    /// \brief Evaluates the constraint jacobian against an input matrix mat.
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
      const Eigen::MatrixBase<OutputMatrix> & res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      derived().jacobianMatrixProductImpl(
        model, data, cdata, mat.derived(), res.const_cast_derived(), aot);
    }

    /// \brief Evaluates the transpose of the constraint jacobian against an input matrix mat.
    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProduct(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      return derived().jacobianTransposeMatrixProductImpl(model, data, cdata, mat.derived());
    }

    /// \brief Evaluates the transpose of the constraint jacobian against an input matrix mat.
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
      const Eigen::MatrixBase<OutputMatrix> & res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      derived().jacobianTransposeMatrixProductImpl(
        model, data, cdata, mat.derived(), res.const_cast_derived(), aot);
    }

    /// \brief Map the constraint forces (aka constraint Lagrange multipliers) to joint space (e.g.,
    /// joint forces, joint torque vector).
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] constraint_forces Input constraint forces (Lagrange multipliers) associated with
    /// the constraint.
    /// \param[out] joint_forces Output joint forces associated with each joint of the model.
    /// \param[out] joint_torques Output joint torques associated with the model.
    /// \param[in] reference_frame Input reference frame in which the forces are expressed.
    ///
    /// \note The results will be added to the joint_forces and joint_torques ouput argument.
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForceLike,
      typename ForceAllocator,
      typename JointTorquesLike,
      ReferenceFrame rf>
    void mapConstraintForceToJointSpace(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques,
      ReferenceFrameTag<rf> reference_frame) const
    {
      derived().mapConstraintForceToJointSpaceImpl(
        model, data, cdata, constraint_forces, joint_forces, joint_torques.const_cast_derived(),
        reference_frame);
    }

    /// \brief Map the joint space quantities (e.g.,
    /// joint motions, joint motion vector) to the constraint motions.
    ///
    /// \param[in] model The model of the rigid body system.
    /// \param[in] data The data associated with model.
    /// \param[in] cdata The constraint data associated with the constraint model.
    /// \param[in] joint_motions Input joint motions associated with the model.
    /// \param[in] joint_generalized_velocity Input joint motions associated with the model.
    /// \param[out] constraint_motions Output constraint motions.
    /// \param[in] reference_frame Input reference frame in which the joint motions are expressed.
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename JointMotionsLike,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointSpaceToConstraintMotion(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<JointMotionsLike> & joint_generalized_velocity,
      const Eigen::MatrixBase<VectorLike> & constraint_motions,
      ReferenceFrameTag<rf> reference_frame) const
    {
      derived().mapJointSpaceToConstraintMotionImpl(
        model, data, cdata, joint_motions, joint_generalized_velocity, constraint_motions,
        reference_frame);
    }

    /// \brief Append to data the apparent inertia due to the constraint.
    template<
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertias(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      derived().appendCouplingConstraintInertiasImpl(
        model, data, cdata, diagonal_constraint_inertia.derived(), reference_frame);
    }

    /// \brief Returns the colwise sparsity associated with a given row
    const BooleanVector & getRowSparsityPattern(const Eigen::Index row_id) const
    {
      return derived().getRowSparsityPatternImpl(row_id);
    }

    /// \brief Returns the colwise sparsity associated with a given row of the active set of
    /// the constraints.
    /// \note If constraints are dynamic (e.g. joint limits), this vector is computed when
    /// calling the calc method.
    const BooleanVector & getActiveRowSparsityPattern(
      const ConstraintData & constraint_data, const Eigen::Index row_id) const
    {
      if constexpr (traits<Derived>::constant_size)
      {
        return getRowSparsityPattern(row_id);
      }
      else
      {
        return derived().getActiveRowSparsityPatternImpl(constraint_data, row_id);
      }
    }

    /// \brief Returns the vector of the activable indexes associated with a given row
    const EigenIndexVector & getActivableRowIndexes(const Eigen::DenseIndex row_id) const
    {
      return derived().getActivableRowIndexesImpl(row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    /// \note If constraints are dynamic (e.g. joint limits), this vector is computed when
    /// calling the calc method.
    const EigenIndexVector & getActiveRowIndexes(
      const ConstraintData & constraint_data, const Eigen::DenseIndex row_id) const
    {
      if constexpr (traits<Derived>::constant_size)
      {
        return getActivableRowIndexes(row_id);
      }
      else
      {
        return derived().getActiveRowIndexesImpl(constraint_data, row_id);
      }
    }

    /// \brief Returns the active compliance internally stored in the constraint and corresponding
    /// to the active set contained in cdata
    /// \note If constraints are dynamic (e.g. joint limits), this vector is computed when
    /// calling the calc method.
    ActiveComplianceVectorTypeConstRef
    getActiveCompliance(const ConstraintData & constraint_data) const
    {
      if constexpr (traits<Derived>::constant_size)
      {
        return compliance();
      }
      else
      {
        return derived().getActivecomplianceImpl(constraint_data);
      }
    }

    /// \brief Returns the active compliance internally stored in the constraint and corresponding
    /// to the active set contained in cdata
    /// \note If constraints are dynamic (e.g. joint limits), this vector is computed when
    /// calling the calc method.
    ActiveComplianceVectorTypeRef getActiveCompliance(ConstraintData & constraint_data) const
    {
      if constexpr (traits<Derived>::constant_size)
      {
        return compliance();
      }
      else
      {
        return derived().getActivecomplianceImpl(constraint_data);
      }
    }

    /// \brief Returns the compliance internally stored in the constraint model.
    ComplianceVectorTypeConstRef compliance() const
    {
      return derived().compliance_impl();
    }

    /// \brief Returns the compliance internally stored in the constraint model.
    ComplianceVectorTypeRef compliance()
    {
      return derived().compliance_impl();
    }

    // CHOICE: right now we use the scalar Baumgarte
    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersConstRef baumgarte_corrector_vector_parameters() const
    // {
    //   return derived().baumgarte_corrector_vector_parameters_impl();
    // }

    // /// \brief Returns the Baumgarte vector parameters internally stored in the constraint model
    // BaumgarteCorrectorVectorParametersRef baumgarte_corrector_vector_parameters()
    // {
    //   return derived().baumgarte_corrector_vector_parameters_impl();
    // }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    const BaumgarteCorrectorParameters & baumgarte_corrector_parameters() const
    {
      return derived().baumgarte_corrector_parameters_impl();
    }

    /// \brief Returns the Baumgarte parameters internally stored in the constraint model
    BaumgarteCorrectorParameters & baumgarte_corrector_parameters()
    {
      return derived().baumgarte_corrector_parameters_impl();
    }

    /// \brief Returns an instance of the associated constraint set operator.
    ConstraintSet set() const
    {
      return derived().setImpl();
    }

    // Attributes common to all constraints

    /// \brief Name of the constraint
    std::string name;

  protected:
    /// \brief Constructor from model.
    /// Protected so that ConstraintModelBase cannot be constructed.
    template<int Options, template<typename, int> class JointCollectionTpl>
    explicit ConstraintModelBase(const ModelTpl<Scalar, Options, JointCollectionTpl> & /*model*/)
    {
    }

    /// \brief Default constructor
    /// Protected so that ConstraintModelBase cannot be constructed.
    ConstraintModelBase()
    {
    }
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_constraint_model_base_hpp__
