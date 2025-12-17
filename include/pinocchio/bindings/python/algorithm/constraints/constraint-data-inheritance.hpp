//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
#define __pinocchio_python_algorithm_constraints_data_inheritance_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<class T>
    struct BinaryKinematicsConstraintDataBasePythonVisitor
    : public bp::def_visitor<BinaryKinematicsConstraintDataBasePythonVisitor<T>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<>(bp::arg("self"), "Default constructor."))
          .def(bp::init<const typename T::ConstraintModel &>(
            bp::args("self", "constraint_model"), "From model constructor."))
          .PINOCCHIO_ADD_PROPERTY(T, constraint_force, "Resulting force.")
          .PINOCCHIO_ADD_PROPERTY(T, oMc1, "Placement of the constraint frame 1 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(T, oMc2, "Placement of the constraint frame 2 wrt WORLD.")
          .PINOCCHIO_ADD_PROPERTY(T, c1Mc2, "Placement of the constraint frame 2 wrt frame 1.")
          .PINOCCHIO_ADD_PROPERTY(T, constraint_position_error, "Constraint position error.")
          .PINOCCHIO_ADD_PROPERTY(T, constraint_velocity_error, "Constraint velocity error.")
          .PINOCCHIO_ADD_PROPERTY(
            T, constraint_acceleration_error, "Constraint acceleration error.")
          .PINOCCHIO_ADD_PROPERTY(
            T, constraint_acceleration_biais_term, "Constraint acceleration term.");
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
