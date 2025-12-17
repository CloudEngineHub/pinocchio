//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
#define __pinocchio_python_algorithm_constraints_model_inheritance_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // BinaryKinematicsConstraintModelBasePythonVisitor
    template<class T>
    struct BinaryKinematicsConstraintModelBasePythonVisitor
    : public bp::def_visitor<BinaryKinematicsConstraintModelBasePythonVisitor<T>>
    {
      typedef typename T::Scalar Scalar;
      typedef typename T::ConstraintSet ConstraintSet;
      typedef typename T::ConstraintData ConstraintData;
      typedef typename T::MatrixSize6 MatrixSize6;
      typedef typename T::SE3 SE3;
      typedef context::Model Model;
      typedef context::Data Data;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(
            bp::init<const Model &, JointIndex, const SE3 &, JointIndex, const SE3 &>(
              (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint1_placement"),
               bp::arg("joint2_id"), bp::arg("joint2_placement")),
              "Contructor from given joint index and placement for the two joints "
              "implied in the constraint."))
          .def(
            bp::init<const Model &, JointIndex>(
              (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id")),
              "Contructor from given joint index of the first joint "
              "implied in the constraint."))
          .def(
            bp::init<const Model &, JointIndex, const SE3 &>(
              (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"),
               bp::arg("joint1_placement")),
              "Contructor from given joint index and placement of the first joint "
              "implied in the constraint."))
          .def(
            bp::init<const Model &, JointIndex, JointIndex>(
              (bp::arg("self"), bp::arg("model"), bp::arg("joint1_id"), bp::arg("joint2_id")),
              "Contructor from given joint index for the two joints "
              "implied in the constraint."))
          .def(
            "getA1", &getA1, bp::args("self", "constraint_data", "reference_frame"),
            "Returns the constraint projector associated with joint 1. "
            "This matrix transforms a spatial velocity expressed in a reference frame "
            "to the first component of the constraint associated with joint 1.")
          .def(
            "getA2", &getA2, bp::args("self", "constraint_data", "reference_frame"),
            "Returns the constraint projector associated with joint 2. "
            "This matrix transforms a spatial velocity expressed in a reference frame "
            "to the first component of the constraint associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(T, joint1_id, "Index of the first joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(T, joint2_id, "Index of the second joint in the model tree.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint1_placement, "Position of attached point with respect to the frame of joint1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint2_placement, "Position of attached point with respect to the frame of joint2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_offset, "Desired constraint shift at position level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_velocity, "Desired constraint velocity at velocity level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, desired_constraint_acceleration,
            "Desired constraint velocity at acceleration level.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_joint1_sparsity, "Colwise sparsity pattern associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_joint2_sparsity, "Colwise sparsity pattern associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint1_span_indexes, "Jointwise span indexes associated with joint 1.")
          .PINOCCHIO_ADD_PROPERTY(
            T, joint2_span_indexes, "Jointwise span indexes associated with joint 2.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_sparsity, "Sparsity pattern associated to the constraint.")
          .PINOCCHIO_ADD_PROPERTY(
            T, colwise_span_indexes, "Indexes of the columns spanned by the constraints.");
      }

      static MatrixSize6
      getA1(const T & self, const ConstraintData & constraint_data, ReferenceFrame rf)
      {
        MatrixSize6 res;
        switch (rf)
        {
        case WORLD:
          res = self.getA1(constraint_data, WorldFrameTag());
        case LOCAL:
          res = self.getA1(constraint_data, LocalFrameTag());
        case LOCAL_WORLD_ALIGNED:
          res = self.getA1(constraint_data, LocalWorldAlignedFrameTag());
        }
        return res;
      }

      static MatrixSize6
      getA2(const T & self, const ConstraintData & constraint_data, ReferenceFrame rf)
      {
        MatrixSize6 res;
        switch (rf)
        {
        case WORLD:
          res = self.getA2(constraint_data, WorldFrameTag());
        case LOCAL:
          res = self.getA2(constraint_data, LocalFrameTag());
        case LOCAL_WORLD_ALIGNED:
          res = self.getA2(constraint_data, LocalWorldAlignedFrameTag());
        }
        return res;
      }
    }; // BinaryKinematicsConstraintModelBasePythonVisitor

    // PointConstraintModelBasePythonVisitor
    template<class T>
    struct PointConstraintModelBasePythonVisitor
    : public bp::def_visitor<PointConstraintModelBasePythonVisitor<T>>
    {
      typedef typename T::Scalar Scalar;
      typedef typename T::ConstraintSet ConstraintSet;
      typedef typename T::ConstraintData ConstraintData;
      typedef typename T::MatrixSize6 MatrixSize6;
      typedef context::Model Model;
      typedef context::Data Data;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(
          "computeConstraintSpatialInertia", &computeConstraintSpatialInertia,
          bp::args("self", "placement", "diagonal_constraint_inertia"),
          "This function computes the spatial inertia associated with the constraint.");
        // computeConstraintInertias is not exposed as it is designed for Eigen Blocks
      }

      static context::Matrix6s computeConstraintSpatialInertia(
        const T & self,
        const context::SE3 & placement,
        const context::Vector3s & diagonal_constraint_inertia)
      {
        return self.computeConstraintSpatialInertia(placement, diagonal_constraint_inertia);
      }
    }; // PointConstraintModelBasePythonVisitor

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
