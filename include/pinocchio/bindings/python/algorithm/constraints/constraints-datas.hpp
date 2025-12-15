//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_datas_hpp__
#define __pinocchio_python_algorithm_constraints_datas_hpp__

#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-data-inheritance.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // generic expose_constraint_data : do nothing special
    template<class T>
    inline bp::class_<T> & expose_constraint_data(bp::class_<T> & cl)
    {
      return cl;
    }

    // specialization for ConstraintDatas
    template<>
    bp::class_<context::PointAnchorConstraintData> &
    expose_constraint_data(bp::class_<context::PointAnchorConstraintData> & cl)
    {
      return cl.def(
        BinaryKinematicsConstraintDataBasePythonVisitor<context::PointAnchorConstraintData>());
    }

    template<>
    bp::class_<context::PointContactData> &
    expose_constraint_data(bp::class_<context::PointContactData> & cl)
    {
      return cl.def(BinaryKinematicsConstraintDataBasePythonVisitor<context::PointContactData>());
    }

    template<>
    bp::class_<context::FrameAnchorConstraintData> &
    expose_constraint_data(bp::class_<context::FrameAnchorConstraintData> & cl)
    {
      return cl.def(
        BinaryKinematicsConstraintDataBasePythonVisitor<context::FrameAnchorConstraintData>());
    }

    template<>
    bp::class_<context::JointFrictionConstraintData> &
    expose_constraint_data(bp::class_<context::JointFrictionConstraintData> & cl)
    {
      return cl.def(
        bp::init<const typename context::JointFrictionConstraintData::ConstraintModel &>(
          bp::args("self", "constraint_model"), "From model constructor."));
    }

    template<>
    bp::class_<context::JointLimitConstraintData> &
    expose_constraint_data(bp::class_<context::JointLimitConstraintData> & cl)
    {
      typedef context::JointLimitConstraintData Self;
      return cl
        .def(bp::init<const typename Self::ConstraintModel &>(
          bp::args("self", "constraint_model"), "From model constructor."))
        .PINOCCHIO_ADD_PROPERTY(
          Self, active_idx_in_activable, "List of idx in all activable that are active.")
        .PINOCCHIO_ADD_PROPERTY(
          Self, active_idx_in_selected, "List of idx in selected joints for each active.")
        .PINOCCHIO_ADD_PROPERTY(
          Self, active_idx_qs_reduce, "List of idx in [0,nqreduce] for each active.")
        .PINOCCHIO_ADD_PROPERTY(
          Self, lower_residual_size, "Number of lower among active constraints.")
        .PINOCCHIO_ADD_PROPERTY(
          Self, activable_constraint_residual, "Activable constraint residual.")
        .PINOCCHIO_ADD_PROPERTY(Self, rowise_tangent_map, "Rowise tangent map.")
        .add_property(
          "constraint_residual",
          bp::make_function(
            +[](const context::JointLimitConstraintData & self)
              -> Eigen::Ref<context::JointLimitConstraintData::VectorXs> {
              return Eigen::Ref<context::JointLimitConstraintData::VectorXs>(
                self.constraint_residual);
            },
            bp::with_custodian_and_ward_postcall<0, 1>()),
          "");
      // CompactTangentMap and constraint_residual_storage are not exposed
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_datas_hpp__
