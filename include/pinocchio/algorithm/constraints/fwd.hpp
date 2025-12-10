//
// Copyright (c) 2022-2024 INRIA
//

#ifndef __pinocchio_algorithm_constraints_fwd_hpp__
#define __pinocchio_algorithm_constraints_fwd_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include <boost/variant.hpp>

namespace pinocchio
{
  // Constraints
  template<typename Scalar, int Options = 0>
  struct RigidConstraintModelTpl;
  template<typename Scalar, int Options = 0>
  struct RigidConstraintDataTpl;

  template<typename Scalar, int Options = 0>
  struct JointFrictionConstraintModelTpl;
  typedef JointFrictionConstraintModelTpl<context::Scalar> JointFrictionConstraintModel;

  template<typename Scalar, int Options = 0>
  struct JointFrictionConstraintDataTpl;
  typedef JointFrictionConstraintDataTpl<context::Scalar> JointFrictionConstraintData;

  template<typename Scalar, int Options = 0>
  struct JointLimitConstraintModelTpl;
  typedef JointLimitConstraintModelTpl<context::Scalar> JointLimitConstraintModel;

  template<typename Scalar, int Options = 0>
  struct JointLimitConstraintDataTpl;
  typedef JointLimitConstraintDataTpl<context::Scalar> JointLimitConstraintData;

  template<typename Scalar, int Options = 0>
  struct PointAnchorConstraintModelTpl;
  typedef PointAnchorConstraintModelTpl<context::Scalar> PointAnchorConstraintModel;
  template<typename Scalar, int Options = 0>
  struct PointAnchorConstraintDataTpl;
  typedef PointAnchorConstraintDataTpl<context::Scalar> PointAnchorConstraintData;

  template<typename Scalar, int Options = 0>
  struct PointContactModelTpl;
  typedef PointContactModelTpl<context::Scalar> PointContactModel;
  template<typename Scalar, int Options = 0>
  struct PointContactDataTpl;
  typedef PointContactDataTpl<context::Scalar> PointContactData;

  template<typename Scalar, int Options = 0>
  struct FrameAnchorConstraintModelTpl;
  typedef FrameAnchorConstraintModelTpl<context::Scalar> FrameAnchorConstraintModel;
  template<typename Scalar, int Options = 0>
  struct FrameAnchorConstraintDataTpl;
  typedef FrameAnchorConstraintDataTpl<context::Scalar> FrameAnchorConstraintData;

  template<typename Scalar, int Options = 0>
  struct ConstraintCollectionDefaultTpl;

  typedef ConstraintCollectionDefaultTpl<context::Scalar, context::Options>
    ConstraintCollectionDefault;

  template<
    typename Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl = ConstraintCollectionDefaultTpl>
  struct ConstraintModelTpl;
  typedef ConstraintModelTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>
    ConstraintModel;

  template<
    typename Scalar,
    int _Options,
    template<typename S, int O> class ConstraintCollectionTpl = ConstraintCollectionDefaultTpl>
  struct ConstraintDataTpl;
  typedef ConstraintDataTpl<context::Scalar, context::Options, ConstraintCollectionDefaultTpl>
    ConstraintData;

  // Sets
  template<typename Scalar, int Options = 0>
  struct BoxSetTpl;
  typedef BoxSetTpl<context::Scalar> BoxSet;

  // Cone sets
  template<typename Scalar, int Options = 0>
  struct FullSpaceConeTpl;
  typedef FullSpaceConeTpl<context::Scalar> FullSpaceCone;

  template<typename Scalar, int Options = 0>
  struct ZeroConeTpl;
  typedef ZeroConeTpl<context::Scalar> ZeroCone;

  template<typename Scalar>
  struct CoulombFrictionConeTpl;
  typedef CoulombFrictionConeTpl<context::Scalar> CoulombFrictionCone;

  template<typename Scalar>
  struct DualCoulombFrictionConeTpl;
  typedef DualCoulombFrictionConeTpl<context::Scalar> DualCoulombFrictionCone;

  template<typename Scalar>
  struct NonNegativeOrthantConeTpl;
  typedef NonNegativeOrthantConeTpl<context::Scalar> NonNegativeOrthantCone;
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_fwd_hpp__
