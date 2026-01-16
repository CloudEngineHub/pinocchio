//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_point_contact_constraint_hpp__
#define __pinocchio_algorithm_constraints_point_contact_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/sets/coulomb-friction-cone.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, PointContactConstraintModelTpl<Scalar, Options>>
  {
    typedef PointContactConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<PointContactConstraintModelTpl<_Scalar, _Options>>
  : traits<PointConstraintModelBase<PointContactConstraintModelTpl<_Scalar, _Options>>>
  {
    typedef _Scalar Scalar;

    static constexpr int Options = _Options;
    static constexpr int Size = 3;

    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef PointContactConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef PointContactConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;
    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    // --------------------------------------------------------------
    // Traits for the algorithmic methods on current state
    // --------------------------------------------------------------
    // Elementary types
    typedef Eigen::Matrix<Scalar, Size, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Eigen::Matrix<Scalar, Size, 1, Options> VectorConstraintSize;

    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;

    // -------------------------------
    // Traits for holded Data
    // -------------------------------
    typedef Eigen::Matrix<Scalar, Size, 1, Options> ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;
  };

  template<typename _Scalar, int _Options>
  struct traits<PointContactConstraintDataTpl<_Scalar, _Options>>
  : traits<PointContactConstraintModelTpl<_Scalar, _Options>>
  {
  };

  ///
  /// \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PointContactConstraintModelTpl
  : PointConstraintModelBase<PointContactConstraintModelTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    typedef PointConstraintModelBase<PointContactConstraintModelTpl> Base;
    typedef ConstraintModelBase<PointContactConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct PointContactConstraintModelTpl;

    typedef PointContactConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;

    using typename Base::SE3;

    using RootBase::classname;

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

    ///
    /// \brief Default constructor
    ///
    PointContactConstraintModelTpl()
    : Base()
    {
    }

    ///
    /// \brief Contructor from joint indexes and placements.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// \param[in] joint2_placement Placement of the constraint w.r.t the frame of joint2.
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointContactConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : Base(model, joint1_id, joint1_placement, joint2_id, joint2_placement)
    {
    }

    ///
    /// \brief Contructor from joint1_id and placement.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointContactConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : Base(model, joint1_id, joint1_placement)
    {
    }

    ///
    /// \brief Contructor from joint ids.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointContactConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model, joint1_id, joint2_id)
    {
    }

    ///
    /// \brief Contructor from joint1_id.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    ///
    /// \remarks The second joint id (joint2_id) is set to be 0 (corresponding to the index of the
    /// universe).
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointContactConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : Base(model, joint1_id)
    {
    }

    // Operators ---------------------

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, PointContactConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, PointContactConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);
      return res;
    }

    ///
    /// \brief Comparison operator
    ///
    /// \param[in] other Other PointContactConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const PointContactConstraintModelTpl & other) const
    {
      return base() == other.base();
    }

    ///
    /// \brief Oposite of the comparison operator.
    ///
    /// \param[in] other Other PointContactConstraintModelTpl to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    bool operator!=(const PointContactConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    /// Specialized accessors --------

    /// \brief Get the friction coefficient of this contact constraint.
    Scalar getFriction() const
    {
      return m_friction;
    }

    /// \brief Set the friction coefficient of this contact constraint.
    void setFriction(Scalar friction)
    {
      PINOCCHIO_THROW_IF(
        check_expression_if_real<Scalar>(friction < 0), std::runtime_error,
        "friction must be > 0 for contact constraints.");
      m_friction = friction;
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    // General -----------------------

    /// \copydoc RootBase::classname
    static std::string classnameImpl()
    {
      return std::string("PointContactConstraintModel");
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

    // Methods for algorithms --------

    /// \copydoc RootBase::set
    ConstraintSet setImpl(const ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return ConstraintSet(m_friction);
    }

  protected:
    Scalar m_friction = Scalar(0.5);

  }; // struct PointContactConstraintModelTpl<_Scalar,_Options>

  ///
  /// \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PointContactConstraintDataTpl
  : PointConstraintDataBase<PointContactConstraintDataTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    typedef PointContactConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef PointContactConstraintDataTpl ConstraintData;
    typedef PointConstraintDataBase<PointContactConstraintDataTpl> Base;

    using typename Base::SE3;

    using Base::classname;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // CRTP related ------------------

    /// \brief Cast to base class
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    /// \brief Const cast to base class
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    // Constructors ------------------

    /// \brief Default constructor
    PointContactConstraintDataTpl()
    {
    }

    /// \brief Constructor from a constraint_model
    explicit PointContactConstraintDataTpl(const ConstraintModel & constraint_model)
    : Base(constraint_model)
    {
    }

    // Operators ---------------------

    /// \brief Comparison operator
    bool operator==(const PointContactConstraintDataTpl & other) const
    {
      return base() == other.base();
    }

    /// \brief Comparison operator
    bool operator!=(const PointContactConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc Base::classname
    static std::string classnameImpl()
    {
      return std::string("PointContactConstraintData");
    }

    /// \copydoc Base::shortname
    std::string shortnameImpl() const
    {
      return classname();
    }
  }; // struct PointContactConstraintDataTpl<_Scalar,_Options>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_point_contact_constraint_hpp__
