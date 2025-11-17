//
// Copyright (c) 2019-2024 INRIA CNRS
//

#ifndef __pinocchio_algorithm_constraints_frame_anchor_constraint_hpp__
#define __pinocchio_algorithm_constraints_frame_anchor_constraint_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/full-space-cone.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-data-base.hpp"

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
      Options = _Options
    };

    typedef FrameAnchorConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef FrameAnchorConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef FullSpaceConeTpl<Scalar, Options> ConstraintSet;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef Eigen::Matrix<Scalar, 6, 1, Options> Vector6;
    typedef Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options> JacobianMatrixType;
    typedef Vector6 VectorConstraintSize;

    typedef Vector6 ComplianceVectorType;
    typedef ComplianceVectorType & ComplianceVectorTypeRef;
    typedef const ComplianceVectorType & ComplianceVectorTypeConstRef;

    typedef ComplianceVectorTypeRef ActiveComplianceVectorTypeRef;
    typedef ComplianceVectorTypeConstRef ActiveComplianceVectorTypeConstRef;

    typedef Vector6 BaumgarteVectorType;
    typedef BaumgarteCorrectorVectorParametersTpl<BaumgarteVectorType>
      BaumgarteCorrectorVectorParameters;
    typedef BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersRef;
    typedef const BaumgarteCorrectorVectorParameters & BaumgarteCorrectorVectorParametersConstRef;
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

    typedef FrameConstraintModelBase<FrameAnchorConstraintModelTpl> Base;
    typedef ConstraintModelBase<FrameAnchorConstraintModelTpl> RootBase;

    template<typename NewScalar, int NewOptions>
    friend struct FrameAnchorConstraintModelTpl;

    typedef FrameAnchorConstraintDataTpl<Scalar, Options> ConstraintData;
    typedef FullSpaceConeTpl<Scalar, Options> ConstraintSet;

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
    FrameAnchorConstraintModelTpl()
    : Base()
    {
    }

    ///
    ///  \brief Contructor with from a given type, joint indexes and placements.
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
    ///  \brief Contructor with from a given type, joint1_id and placement.
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
    ///  \brief Contructor with from a given type and the joint ids.
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
    ///  \brief Contructor with from a given type and .
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

    ///
    /// \brief Create data storage associated to the constraint
    ///
    ConstraintData createData() const
    {
      return ConstraintData(*this);
    }

    /// \brief Cast operator
    template<typename NewScalar>
    typename CastType<NewScalar, FrameAnchorConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, FrameAnchorConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);
      res.m_set = m_set.template cast<NewScalar>();
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
      return base() == other.base() && m_set == other.m_set;
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
      return std::string("FrameAnchorConstraintModel");
    }
    std::string shortname() const
    {
      return classname();
    }

  protected:
    ConstraintSet m_set = ConstraintSet(6);

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

    using typename Base::SE3;

    /// \brief Default constructor
    FrameAnchorConstraintDataTpl()
    {
    }

    explicit FrameAnchorConstraintDataTpl(const ConstraintModel & constraint_model)
    : Base(constraint_model)
    {
    }

    bool operator==(const FrameAnchorConstraintDataTpl & other) const
    {
      return base() == other.base();
    }

    bool operator!=(const FrameAnchorConstraintDataTpl & other) const
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
      return std::string("FrameAnchorConstraintData");
    }
    std::string shortname() const
    {
      return classname();
    }
  }; // struct FrameAnchorConstraintDataTpl<_Scalar,_Options>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_frame_anchor_constraint_hpp__
