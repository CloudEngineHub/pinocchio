//
// Copyright (c) 2025 INRIA
//

#pragma once

#include "pinocchio/algorithm/constraints.hpp"
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
        .def(
          bp::init<const typename Self::ConstraintModel &>(
            bp::args("self", "constraint_model"), "From model constructor."))
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
