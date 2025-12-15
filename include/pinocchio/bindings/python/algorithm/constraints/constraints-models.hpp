//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_models_hpp__
#define __pinocchio_python_algorithm_constraints_models_hpp__

#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-model-inheritance.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // generic expose_constraint_model : do nothing special
    template<class T>
    bp::class_<T> & expose_constraint_model(bp::class_<T> & cl)
    {
      return cl;
    }

    // specialization for ConstraintModels
    template<>
    bp::class_<context::PointAnchorConstraintModel> &
    expose_constraint_model(bp::class_<context::PointAnchorConstraintModel> & cl)
    {
      return cl
        .def(
          BinaryKinematicsConstraintModelBasePythonVisitor<context::PointAnchorConstraintModel>())
        .def(PointConstraintModelBasePythonVisitor<context::PointAnchorConstraintModel>());
    }

    template<>
    bp::class_<context::PointContactModel> &
    expose_constraint_model(bp::class_<context::PointContactModel> & cl)
    {
      typedef context::PointContactModel Self;
      return cl.def(BinaryKinematicsConstraintModelBasePythonVisitor<context::PointContactModel>())
        .def(PointConstraintModelBasePythonVisitor<context::PointContactModel>())
        .def("getFriction", &Self::getFriction, "Get coulomb friction coefficient.")
        .def("setFriction", &Self::setFriction, "Set coulomb friction coefficient.");
    }

    template<>
    bp::class_<context::FrameAnchorConstraintModel> &
    expose_constraint_model(bp::class_<context::FrameAnchorConstraintModel> & cl)
    {
      return cl.def(
        BinaryKinematicsConstraintModelBasePythonVisitor<context::FrameAnchorConstraintModel>());
    }

    template<>
    bp::class_<context::JointFrictionConstraintModel> &
    expose_constraint_model(bp::class_<context::JointFrictionConstraintModel> & cl)
    {
      typedef typename context::JointFrictionConstraintModel::JointIndexVector JointIndexVector;
      typedef context::JointFrictionConstraintModel Self;
      typedef context::VectorXs VectorXs;
      cl.def(bp::init<const context::Model &, const JointIndexVector &>(
               (bp::arg("self"), bp::arg("model"), bp::arg("active_joints")),
               "Contructor from given joint index vector "
               "implied in the constraint."))
        .def(
          "getActiveJoints", &Self::getActiveJoints,
          bp::return_value_policy<bp::copy_const_reference>())
        .def(
          "getActiveDofs", &Self::getActiveDofs,
          bp::return_value_policy<bp::copy_const_reference>())
        .def(
          "getFrictionLowerLimit", &Self::getFrictionLowerLimit,
          bp::return_value_policy<bp::copy_const_reference>(), "Get friction lower limit.")
        .def(
          "setFrictionLowerLimit", bp::make_function(+[](Self & self, const VectorXs & lb) {
            self.setFrictionLowerLimit(lb);
          }),
          "Set friction lower limit.")
        .def(
          "getFrictionUpperLimit", &Self::getFrictionUpperLimit,
          bp::return_value_policy<bp::copy_const_reference>(), "Get friction upper limit.")
        .def(
          "setFrictionUpperLimit", bp::make_function(+[](Self & self, const VectorXs & ub) {
            self.setFrictionUpperLimit(ub);
          }),
          "Set friction upper limit.");
      return cl;
    }

    template<>
    bp::class_<context::JointLimitConstraintModel> &
    expose_constraint_model(bp::class_<context::JointLimitConstraintModel> & cl)
    {
      typedef typename context::JointLimitConstraintModel::JointIndexVector JointIndexVector;
      typedef typename context::JointLimitConstraintModel Self;
      typedef typename context::JointLimitConstraintData ConstraintData;
      cl.def(bp::init<const context::Model &, const JointIndexVector &>(
               (bp::arg("self"), bp::arg("model"), bp::arg("activable_joints")),
               "Contructor from given joint index vector "
               "implied in the constraint."))
        .def(bp::init<
             const context::Model &, const JointIndexVector &, const context::VectorXs &,
             const context::VectorXs &>(
          (bp::arg("self"), bp::arg("model"), bp::arg("activable_joints"), bp::arg("lb"),
           bp::arg("ub")),
          "Contructor from given joint index vector "
          "implied in the constraint."))
        .def(bp::init<
             const context::Model &, const JointIndexVector &, const context::VectorXs &,
             const context::VectorXs &, const context::VectorXs &>(
          (bp::arg("self"), bp::arg("model"), bp::arg("activable_joints"), bp::arg("lb"),
           bp::arg("ub"), bp::arg("margin")),
          "Contructor from given joint index vector "
          "implied in the constraint."))
        .def(
          "getSelectedJoints", &Self::getSelectedJoints,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Joints for which there is at least one position limit.")
        .def("getNqReduce", &Self::getNqReduce, "Sum of nq of activable joints.")
        .def("getMaxOfNvs", &Self::getMaxOfNvs, "Max nv of atomic joints in activable joints.")
        .def(
          "getActivablePositionLimit", &Self::getActivablePositionLimit,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Position limit of the dof of the constraints.")
        .def(
          "getActivablePositionMargin", &Self::getActivablePositionMargin,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Position margin of the dof of the constraints.")
        .def(
          "lowerMaxResidualSize", &Self::lowerMaxResidualSize,
          "Part of maxResidualSize() that are lower bound limits.")
        .def(
          "upperMaxResidualSize", &Self::upperMaxResidualSize,
          "Part of maxResidualSize() that are upper bound limits.")
        .def(
          "lowerResidualSize",
          bp::make_function(+[](const Self & self, const ConstraintData & cdata) -> int {
            return self.lowerResidualSize(cdata);
          }),
          "Give the size of constraint that are lower bound in a given state.")
        .def(
          "upperResidualSize",
          bp::make_function(+[](const Self & self, const ConstraintData & cdata) -> int {
            return self.upperResidualSize(cdata);
          }),
          "Give the size of constraint that are upper bound in a given state.")
        .def(
          "setPositionLimitAndMargin",
          bp::make_function(
            +[](
               Self & self, const context::VectorXs & lb, const context::VectorXs & ub,
               const context::VectorXs & margin) -> void {
              self.setPositionLimitAndMargin(lb, ub, margin);
            }),
          "Set position limit and margin for activable constraints from lower_bound, upper_bound "
          "and margin of size model.nq.");
      return cl;
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_models_hpp__
