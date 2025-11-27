//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_frame_anchor_constraint_hpp__
#define __pinocchio_algorithm_constraints_frame_anchor_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/sets/full-space-cone.hpp"

namespace pinocchio
{

  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, FrameAnchorConstraintModelTpl<Scalar, Options>>
  {
    typedef FrameAnchorConstraintModelTpl<NewScalar, Options> type;
  };

  template<typename _Scalar, int _Options>
  struct traits<FrameAnchorConstraintModelTpl<_Scalar, _Options>>
  : traits<FrameConstraintModelBase<FrameAnchorConstraintModelTpl<_Scalar, _Options>>>
  {
    typedef _Scalar Scalar;

    enum
    {
      Options = _Options,
      Size = 6
    };

    // --------------------------------------------------------------
    // Traits referencing the constraint and associated types
    // --------------------------------------------------------------
    typedef FrameAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrameAnchorConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef FullSpaceConeTpl<Scalar, Options> ConstraintSet;
    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    // --------------------------------------------------------------
    // Traits for the algorithmic methods on current state
    // --------------------------------------------------------------
    // Elementary types
    typedef Eigen::Matrix<Scalar, Size, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Eigen::Matrix<Scalar, Size, 1, Options> VectorConstraintSize;

    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;

    // -------------------------------
    // Traits for holded Data
    // -------------------------------
    typedef Eigen::Matrix<Scalar, Size, 1, Options> ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    // Will be removed
    typedef ComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef ComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;

    // Not needed anymore
    // typedef Vector3 BaumgarteVectorType;
    // typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorType>
    //   BaumgarteCorrectorVectorParameters;
    // typedef BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersRef;
    // typedef const BaumgarteCorrectorVectorParameters &
    // BaumgarteCorrectorVectorParametersConstRef;
  };

  template<typename _Scalar, int _Options>
  struct traits<FrameAnchorConstraintDataTpl<_Scalar, _Options>>
  : traits<FrameAnchorConstraintModelTpl<_Scalar, _Options>>
  {
  };

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct FrameAnchorConstraintModelTpl
  : FrameConstraintModelBase<FrameAnchorConstraintModelTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef FrameAnchorConstraintModelTpl Self;
    typedef FrameConstraintModelBase<FrameAnchorConstraintModelTpl> Base;
    typedef ConstraintModelBase<FrameAnchorConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct FrameAnchorConstraintModelTpl;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::ConstraintSet ConstraintSet;

    using typename Base::SE3;

    using RootBase::classname;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

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

    ///
    /// \brief Default constructor
    ///
    FrameAnchorConstraintModelTpl()
    : Base()
    {
    }

    ///
    /// \brief Contructor from joint indexes and placements.
    ///
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// \param[in] joint2_placement Placement of the constraint w.r.t the frame of joint2.
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrameAnchorConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : Base(model, joint1_id, joint1_placement, joint2_id, joint2_placement)
    {
    }

    ///
    ///  \brief Contructor from joint1_id and placement.
    ///
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrameAnchorConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : Base(model, joint1_id, joint1_placement)
    {
    }

    ///
    ///  \brief Contructor from joint ids.
    ///
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrameAnchorConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model, joint1_id, joint2_id)
    {
    }

    ///
    ///  \brief Contructor with from joint1_id.
    ///
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    ///
    /// \remarks The second joint id (joint2_id) is set to be 0 (corresponding to the index of the
    /// universe).
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    FrameAnchorConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : Base(model, joint1_id)
    {
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, FrameAnchorConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, FrameAnchorConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);
      return res;
    }

    ///
    ///  \brief Comparison operator
    ///
    /// \param[in] other Other FrameAnchorConstraintModelTpl to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const FrameAnchorConstraintModelTpl & other) const
    {
      return base() == other.base();
    }

    ///
    ///  \brief Oposite of the comparison operator.
    ///
    /// \param[in] other Other FrameAnchorConstraintModelTpl to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    bool operator!=(const FrameAnchorConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc RootBase::createData
    ConstraintData createDataImpl() const
    {
      return ConstraintData(*this);
    }

    /// \copydoc RootBase::set
    ConstraintSet setImpl() const
    {
      return ConstraintSet();
    }

    /// \copydoc RootBase::classname
    static std::string classnameImpl()
    {
      return std::string("FrameAnchorConstraintModel");
    }

    /// \copydoc RootBase::shortname
    std::string shortnameImpl() const
    {
      return classname();
    }

  }; // struct FrameAnchorConstraintModelTpl<_Scalar,_Options>

  ///
  ///  \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename _Scalar, int _Options>
  struct FrameAnchorConstraintDataTpl
  : FrameConstraintDataBase<FrameAnchorConstraintDataTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef FrameAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrameAnchorConstraintDataTpl ConstraintData;
    typedef FrameConstraintDataBase<FrameAnchorConstraintDataTpl> Base;
    typedef ConstraintDataBase<FrameAnchorConstraintDataTpl> RootBase;

    using typename Base::SE3;

    using RootBase::classname;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

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

    /// \brief Default constructor
    FrameAnchorConstraintDataTpl()
    {
    }

    /// \brief Constructor from constraint_model
    explicit FrameAnchorConstraintDataTpl(const ConstraintModel & constraint_model)
    : Base(constraint_model)
    {
    }

    /// \brief Comparison operator
    bool operator==(const FrameAnchorConstraintDataTpl & other) const
    {
      return base() == other.base();
    }

    /// \brief Comparison operator
    bool operator!=(const FrameAnchorConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    /// \copydoc RootBase::classname
    static std::string classnameImpl()
    {
      return std::string("FrameAnchorConstraintData");
    }

    /// \copydoc RootBase::shortname
    std::string shortnameImpl() const
    {
      return classname();
    }
  }; // struct FrameAnchorConstraintDataTpl<_Scalar,_Options>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_frame_anchor_constraint_hpp__
