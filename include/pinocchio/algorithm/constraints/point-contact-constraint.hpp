//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_point_contact_constraint_hpp__
#define __pinocchio_algorithm_constraints_point_contact_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/coulomb-friction-cone.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"

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

    enum
    {
      Options = _Options
    };

    typedef PointContactConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef PointContactConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Vector3 VectorConstraintSize;

    typedef Vector3 ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef ComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef ComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;

    typedef Vector3 BaumgarteVectorType;
    typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorType>
      BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersRef;
    typedef const BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersConstRef;
  };

  template<typename _Scalar, int _Options>
  struct traits<PointContactConstraintDataTpl<_Scalar, _Options>>
  : traits<PointContactConstraintModelTpl<_Scalar, _Options>>
  {
  };

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PointContactConstraintModelTpl
  : PointConstraintModelBase<PointContactConstraintModelTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef PointConstraintModelBase<PointContactConstraintModelTpl> Base;
    typedef ConstraintModelBase<PointContactConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct PointContactConstraintModelTpl;

    typedef PointContactConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef CoulombFrictionConeTpl<Scalar> ConstraintSet;

    using typename Base::SE3;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    ///
    ///  \brief Default constructor
    ///
    PointContactConstraintModelTpl()
    : Base()
    {
    }

    ///
    ///  \brief Contructor with from a given type, joint indexes and placements.
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
    ///  \brief Contructor with from a given type, joint1_id and placement.
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
    ///  \brief Contructor with from a given type and the joint ids.
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
    ///  \brief Contructor with from a given type and .
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

    ///
    /// \brief Create data storage associated to the constraint
    ///
    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

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
    ///  \brief Comparison operator
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
    ///  \brief Oposite of the comparison operator.
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

    /// \brief Get the friction coefficient of this contact constraint.
    Scalar getFriction() const
    {
      return m_friction;
    }

    /// \brief Set the friction coefficient of this contact constraint.
    void setFriction(Scalar friction)
    {
      PINOCCHIO_THROW_IF(
        friction < 0, std::runtime_error, "friction must be > 0 for contact constraints.");
      m_friction = friction;
    }

    /// \copydoc Base::set
    ConstraintSet setImpl() const
    {
      return ConstraintSet(m_friction);
    }

    static std::string classname()
    {
      return std::string("PointContactConstraintModel");
    }
    std::string shortname() const
    {
      return classname();
    }

  protected:
    Scalar m_friction = Scalar(0.5);

  }; // struct PointContactConstraintModelTpl<_Scalar,_Options>

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PointContactConstraintDataTpl
  : PointConstraintDataBase<PointContactConstraintDataTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef PointContactConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef PointContactConstraintDataTpl ConstraintData;
    typedef PointConstraintDataBase<PointContactConstraintDataTpl> Base;

    using typename Base::SE3;

    /// \brief Default constructor
    PointContactConstraintDataTpl()
    {
    }

    explicit PointContactConstraintDataTpl(const ConstraintModel & constraint_model)
    : Base(constraint_model)
    {
    }

    bool operator==(const PointContactConstraintDataTpl & other) const
    {
      return base() == other.base();
    }

    bool operator!=(const PointContactConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    Base & base()
    {
      return static_cast<Base &>(*this);
    }
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    static std::string classname()
    {
      return std::string("PointContactConstraintData");
    }
    std::string shortname() const
    {
      return classname();
    }
  }; // struct PointContactConstraintDataTpl<_Scalar,_Options>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_point_contact_constraint_hpp__
