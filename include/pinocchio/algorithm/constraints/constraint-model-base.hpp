//
// Copyright (c) 2023-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraints_constraint_model_base_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/common/model-entity.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

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

  enum struct ConstraintSizeType
  {
    STATIC,   // The size is fixed at compile time. residualSize = maxResidualSize
              // -> no implementation to provide.
    CONSTANT, // The size is fixed at build time. residualSize = maxResidualSize
              // -> maxResidualSizeImpl to implement
    BOUNDED,  // The maxResidualSize is fixed at build time. 0 <= size <= maxResidualSize
              // -> residualSizeImpl and maxResidualSizeImpl to implement
    GENERAL   // The size is not guaranteed to be bounded
              // -> residualSizeImpl to implement
  };

  template<class Derived>
  struct ConstraintModelBase
  : NumericalBase<Derived>
  , ModelEntity<Derived>
  {
    // --------------------------------------------------------------
    // Type defs
    // --------------------------------------------------------------

    typedef ModelEntity<Derived> Base;

    // Retrieving traits --------------------------------------------
    typedef typename traits<Derived>::Scalar Scalar;
    enum
    {
      Options = traits<Derived>::Options,
      Size = traits<Derived>::Size
    };

    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      traits<Derived>::constraint_formulation_level;
    static constexpr ConstraintSizeType constraint_size_type =
      traits<Derived>::constraint_size_type;
    // {STATIC} \subset {CONSTANT} \subset {BOUNDED} \subset {GENERAL}
    static constexpr bool static_size = constraint_size_type == ConstraintSizeType::STATIC;
    static constexpr bool constant_size =
      static_size || (constraint_size_type == ConstraintSizeType::CONSTANT);
    static constexpr bool bounded_size =
      constant_size || (constraint_size_type == ConstraintSizeType::BOUNDED);

    static constexpr bool has_baumgarte_corrector =
      traits<Derived>::has_baumgarte_corrector; // Baumgarte make sense and exist directly for the
                                                // constraint
    static constexpr bool has_compliance_member =
      traits<Derived>::has_compliance_member; // The constraint itself possesses a member
                                              // m_compliance which can be set by the user
    static constexpr bool has_set =
      traits<Derived>::has_set; // The constraint itself defines the set, otherwise must have a
                                // mechanism for set-related methods

    typedef typename traits<Derived>::ConstraintModel ConstraintModel;
    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ConstraintSet ConstraintSet;

    typedef typename traits<Derived>::JacobianMatrixType JacobianMatrixType;
    typedef typename traits<Derived>::VectorConstraintSize VectorConstraintSize;
    // JacobianMatrixProductReturnType / JacobianTransposeMatrixProductReturnType Tpl

    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef typename traits<Derived>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
    typedef typename traits<Derived>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    // Usefull types ------------------------------------------------
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector;
    typedef std::vector<Eigen::Index> EigenIndexVector;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // CRTP related ------------------

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

    // Constructors ------------------

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

    // Operators ---------------------

  public:
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

    // -------------------------------
    // BASE METHODS
    // -------------------------------

    // General -----------------------

    /// \brief Returns the name of the underlying class if this is a variant.
    static std::string classname()
    {
      return Derived::classnameImpl();
    }

    /// \brief Returns the name of the underlying class if this is a variant.
    std::string shortname() const
    {
      return derived().shortnameImpl();
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

    using Base::createData;

    /// \brief Returns a constraint data associated to this constraint model.
    ConstraintData createData() const
    {
      return derived().createDataImpl();
    }

    // Size management ---------------

    /// \brief Returns the maximum size of the constraint.
    int maxResidualSize() const
    {
      if constexpr (static_size)
      {
        return Size;
      }
      else
      {
        return derived().maxResidualSizeImpl();
      }
    }

    // Methods for algorithm ---------

    /// \brief Returns the current size of the constraint, typically after `calc` has been called.
    /// \note If constraints are dynamic (e.g. joint limits), residualSize is computed when
    /// calling the calc method.
    template<typename ConstraintDataDerived>
    int residualSize(const ConstraintDataBase<ConstraintDataDerived> & cdata) const
    {
      if constexpr (constant_size)
      {
        return maxResidualSize();
      }
      else
      {
        return derived().residualSizeImpl(cdata.derived());
      }
    }

    /// \brief Returns the colwise sparsity associated with a given row of the active set of
    /// the constraints.
    /// \note If constraints are dynamic (e.g. joint limits), this vector is computed when
    /// calling the calc method.
    const BooleanVector &
    getRowSparsityPattern(const ConstraintData & cdata, const Eigen::Index row_id) const
    {
      return derived().getRowSparsityPatternImpl(cdata, row_id);
    }

    /// \brief Returns the vector of the active indexes associated with a given row
    /// \note If constraints are dynamic (e.g. joint limits), this vector is computed when
    /// calling the calc method.
    const EigenIndexVector &
    getRowIndexes(const ConstraintData & cdata, const Eigen::Index row_id) const
    {
      return derived().getRowIndexesImpl(cdata, row_id);
    }

    /// \brief Returns an instance of the associated constraint set operator.
    ConstraintSet set() const
    {
      return derived().setImpl();
    }

    /// \brief Fill the compliance of size residualSize relted to the courant state of the
    /// constraint
    template<typename VectorLike>
    void retrieveCompliance(
      const ConstraintData & cdata, const Eigen::MatrixBase<VectorLike> & res) const
    {
      if constexpr (constant_size && has_compliance_member)
      {
        res.const_cast_derived() = compliance();
      }
      else
      {
        derived().retrieveComplianceImpl(cdata, res.const_cast_derived());
      }
    }

    /// \brief Evaluate the constraint values at the current state given by data and store the
    /// results in cdata.
    /// \note data must be populated by results of a `forwardKinematic(model, data, q, v, a)`.
    /// The forward kinematics on q determines the constraint position error, on v the constraint
    /// velocity error, on a the constraint acceleration error.
    /// Typically, a call to `aba` will fill all the necessary fields of data.
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
      ReturnType res = ReturnType::Zero(residualSize(cdata), model.nv);

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

    // Data handling -----------------

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

    /// \brief Set the compliance
    template<typename VectorLike>
    void setCompliance(const Eigen::MatrixBase<VectorLike> & vector)
    {
      derived().setComplianceImpl(vector);
    }

    /// \brief Set baumgarte corrector
    void setBaumgarteCorrectorParameters(
      const BaumgarteCorrectorParameters & baumgarte_corrector_parameters_in)
    {
      derived().setBaumgarteCorrectorParametersImpl(baumgarte_corrector_parameters_in);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    // General -----------------------
    // classnameImpl()
    // shortnameImpl()
    // createDataImpl()

    // Size Management ---------------
    // maxResidualSizeImpl()  // Not needed for STATIC

    // Methods for algorithms --------
    // residualSizeImpl(const cdata)  // Not needed for <= CONSTANT
    // getRowSparsityPatternImpl(const cdata)
    // getRowIndexesImpl(const cdata)
    // setImpl()  // Not needed if has_set=False
    // retrieveComplianceImpl(cdata, ...)
    // calcImpl(const model, const data, cdata)  // The only one mutating cdata
    // jacobianImpl(const model, const data, const cdata, ...)
    // jacobianMatrixProductImpl(const model, const data, const cdata, ...)
    // jacobianTransposeMatrixProductImpl(const model, const data, const cdata, ...)
    // mapConstraintForceToJointSpaceImpl(const model, const data, const cdata, ...)
    // mapJointSpaceToConstraintMotionImpl(const model, const data, const cdata, ...)
    // appendCouplingConstraintInertiasImpl(const model, const data, const cdata, ...)

    // Data handling -----------------
    // compliance_impl()
    // baumgarte_corrector_parameters_impl()
    // setComplianceImpl(...)
    // setBaumgarteCorrectorParametersImpl(...)

    // Default Implementation --------
    /// \copydoc setCompliance
    template<typename VectorLike>
    void setComplianceImpl(const Eigen::MatrixBase<VectorLike> & vector)
    {
      compliance() = vector;
    }

    /// \copydoc setBaumgarteCorrectorParameters
    void setBaumgarteCorrectorParametersImpl(
      const BaumgarteCorrectorParameters & baumgarte_corrector_parameters_in)
    {
      baumgarte_corrector_parameters() = baumgarte_corrector_parameters_in;
    }

    // ------------------------------
    // MEMBERS
    // ------------------------------

    /// \brief Name of the constraint
    std::string name;
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_constraint_model_base_hpp__
