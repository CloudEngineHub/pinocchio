//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_constraint_collection_default_hpp__
#define __pinocchio_algorithm_constraints_constraint_collection_default_hpp__

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include <boost/variant.hpp>

namespace pinocchio
{
  template<typename _Scalar, int _Options>
  struct ConstraintCollectionDefaultTpl
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };

    typedef PointAnchorConstraintModelTpl<Scalar, Options> PointAnchorConstraintModel;
    typedef PointAnchorConstraintDataTpl<Scalar, Options> PointAnchorConstraintData;

    typedef PointContactModelTpl<Scalar, Options> PointContactModel;
    typedef PointContactDataTpl<Scalar, Options> PointContactData;

    typedef JointFrictionConstraintModelTpl<Scalar, Options> JointFrictionConstraintModel;
    typedef JointFrictionConstraintDataTpl<Scalar, Options> JointFrictionConstraintData;

    typedef JointLimitConstraintModelTpl<Scalar, Options> JointLimitConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> JointLimitConstraintData;

    typedef FrameAnchorConstraintModelTpl<Scalar, Options> FrameAnchorConstraintModel;
    typedef FrameAnchorConstraintDataTpl<Scalar, Options> FrameAnchorConstraintData;

    typedef boost::variant<
      boost::blank,
      PointAnchorConstraintModel,
      PointContactModel,
      JointFrictionConstraintModel,
      JointLimitConstraintModel,
      FrameAnchorConstraintModel>
      ConstraintModelVariant;

    typedef boost::variant<
      boost::blank,
      PointAnchorConstraintData,
      PointContactData,
      JointFrictionConstraintData,
      JointLimitConstraintData,
      FrameAnchorConstraintData>
      ConstraintDataVariant;
  }; // struct ConstraintCollectionDefaultTpl

  typedef ConstraintCollectionDefault::ConstraintModelVariant ConstraintModelVariant;
  typedef ConstraintCollectionDefault::ConstraintDataVariant ConstraintDataVariant;

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_constraint_collection_default_hpp__
