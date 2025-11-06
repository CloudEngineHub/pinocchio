//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_point_anchor_constraint_hpp__
#define __pinocchio_algorithm_constraints_point_anchor_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/unbounded-set.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, PointAnchorConstraintModelTpl<Scalar, Options>>
  {
    typedef PointAnchorConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<PointAnchorConstraintModelTpl<_Scalar, _Options>>
  : traits<PointConstraintModelBase<PointAnchorConstraintModelTpl<_Scalar, _Options>>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options
    };

    typedef PointAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef PointAnchorConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef UnboundedSetTpl<Scalar, Options> ConstraintSet;

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
  struct traits<PointAnchorConstraintDataTpl<_Scalar, _Options>>
  : traits<PointAnchorConstraintModelTpl<_Scalar, _Options>>
  {
  };

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PointAnchorConstraintModelTpl
  : PointConstraintModelBase<PointAnchorConstraintModelTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef PointConstraintModelBase<PointAnchorConstraintModelTpl> Base;
    typedef ConstraintModelBase<PointAnchorConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct PointAnchorConstraintModelTpl;

    typedef PointAnchorConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef UnboundedSetTpl<Scalar, Options> ConstraintSet;

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
    PointAnchorConstraintModelTpl()
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
    PointAnchorConstraintModelTpl(
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
    PointAnchorConstraintModelTpl(
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
    PointAnchorConstraintModelTpl(
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
    PointAnchorConstraintModelTpl(
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
    typename CastType<NewScalar, PointAnchorConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, PointAnchorConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);
      res.m_set = m_set.template cast<NewScalar>();
      return res;
    }

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other PointAnchorConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const PointAnchorConstraintModelTpl & other) const
    {
      return base() == other.base() && m_set == other.m_set;
    }

    ///
    ///  \brief Oposite of the comparison operator.
    ///
    /// \param[in] other Other PointAnchorConstraintModelTpl to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    bool operator!=(const PointAnchorConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    const ConstraintSet & set() const
    {
      return m_set;
    }

    ConstraintSet & set()
    {
      return m_set;
    }

    static std::string classname()
    {
      return std::string("PointAnchorConstraintModel");
    }
    std::string shortname() const
    {
      return classname();
    }

  protected:
    ConstraintSet m_set = ConstraintSet(3);

  }; // struct PointAnchorConstraintModelTpl<_Scalar,_Options>

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct PointAnchorConstraintDataTpl
  : PointConstraintDataBase<PointAnchorConstraintDataTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef PointAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef PointAnchorConstraintDataTpl ConstraintData;
    typedef PointConstraintDataBase<PointAnchorConstraintDataTpl> Base;

    using typename Base::SE3;

    /// \brief Default constructor
    PointAnchorConstraintDataTpl()
    {
    }

    explicit PointAnchorConstraintDataTpl(const ConstraintModel & constraint_model)
    : Base(constraint_model)
    {
    }

    bool operator==(const PointAnchorConstraintDataTpl & other) const
    {
      return base() == other.base();
    }

    bool operator!=(const PointAnchorConstraintDataTpl & other) const
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
      return std::string("PointAnchorConstraintData");
    }
    std::string shortname() const
    {
      return classname();
    }
  }; // struct PointAnchorConstraintDataTpl<_Scalar,_Options>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_point_anchor_constraint_hpp__
