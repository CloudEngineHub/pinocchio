//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
#define __pinocchio_python_algorithm_constraints_model_base_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"
#include "pinocchio/bindings/python/utils/eigen.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/baumgarte-corrector-parameters.hpp"

namespace pinocchio
{
  namespace python
  {

    namespace bp = boost::python;

    template<class ConstraintModelDerived>
    struct ConstraintModelBasePythonVisitor
    : public bp::def_visitor<ConstraintModelBasePythonVisitor<ConstraintModelDerived>>
    {
      typedef ConstraintModelDerived Self;
      typedef typename Self::Scalar Scalar;
      typedef typename Self::ConstraintSet ConstraintSet;
      typedef typename Self::ConstraintData ConstraintData;
      typedef context::Model Model;
      typedef context::Data Data;
      typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
      typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;
      typedef typename traits<Self>::JacobianMatrixType JacobianMatrixType;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.PINOCCHIO_ADD_PROPERTY(Self, name, "Name of the constraint.")
          .def("classname", &Self::classname)
          .staticmethod("classname")
          .def("shortname", &Self::shortname, "Short name of the class.")
          .def(
            "createData", &Self::createData, "Create a Data object for the given constraint model.")
          .def("set", &Self::set, "Constraint set.")
          .add_property(
            "compliance",
            bp::make_function( //
              +[](const Self & self) -> context::VectorXs { return self.compliance(); }),
            bp::make_function( //
              +[](Self & self, const context::VectorXs & new_vector) {
                self.compliance() = new_vector;
              }),
            "Compliance of the constraint.")
          .def(
            "maxSize", +[](const Self & self) -> int { return self.maxSize(); }, bp::arg("self"),
            "Constraint max size.")
          .def(
            "activeSize",
            +[](const Self & self, const ConstraintData & cdata) -> int {
              return self.activeSize(cdata);
            },
            bp::args("self", "constraint_data"), "Constraint state size.")
          .def(
            "calc", &calc, bp::args("self", "model", "data", "constraint_data"),
            "Evaluate the constraint values at the current state given by data and store the "
            "results.")
          .def(
            "jacobian",
            (JacobianMatrixType(Self::*)(const Model &, const Data &, ConstraintData &)
               const)&Self::jacobian,
            bp::args("self", "model", "data", "constraint_data"),
            "Compute the constraint jacobian.")
          .def(
            "jacobianMatrixProduct", &jacobianMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"),
            "Forward chain rule: return product between the jacobian and a matrix.")
          .def(
            "jacobianTransposeMatrixProduct", &jacobianTransposeMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"),
            "Backward chain rule: return product between the jacobian transpose and a matrix.")
          .def(
            "getRowSparsityPattern", &Self::getRowSparsityPattern,
            bp::args("self", "constraint_data", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Active colwise sparsity associated with a given row.")
          .def(
            "getRowIndexes", &Self::getRowIndexes, bp::args("self", "constraint_data", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Vector of the active indexes associated with a given row.")
          .def(
            "getActiveCompliance",
            bp::make_function(
              +[](const Self & self, const ConstraintData & cdata) -> context::VectorXs {
                return self.getActiveCompliance(cdata);
              }),
            "Vector of the active compliance internally stored in the constraint.")
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
        if (::pinocchio::traits<ConstraintModelDerived>::has_baumgarte_corrector)
        {
          typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
          BaumgarteCorrectorParametersPythonVisitor<BaumgarteCorrectorParameters>::expose();

          cl.add_property(
            "baumgarte_corrector_parameters",
            bp::make_function( //
              +[](Self & self) -> BaumgarteCorrectorParameters & {
                return self.baumgarte_corrector_parameters();
              },
              bp::return_internal_reference<>()),
            bp::make_function( //
              +[](Self & self, const BaumgarteCorrectorParameters & copy) {
                self.baumgarte_corrector_parameters() = copy;
              },
              bp::return_internal_reference<>()),
            "Baumgarte parameters associated with the constraint.");
        }
      }

      static void
      resize(Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        self.resize(model, data, constraint_data);
      }

      static void calc(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        self.calc(model, data, constraint_data);
      }

      static context::MatrixXs jacobianMatrixProduct(
        const Self & self,
        const Model & model,
        const Data & data,
        const ConstraintData & constraint_data,
        const context::MatrixXs & matrix)
      {
        context::MatrixXs res =
          context::MatrixXs::Zero(self.activeSize(constraint_data), matrix.cols());
        self.jacobianMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }

      static context::MatrixXs jacobianTransposeMatrixProduct(
        const Self & self,
        const Model & model,
        const Data & data,
        const ConstraintData & constraint_data,
        const context::MatrixXs & matrix)
      {
        context::MatrixXs res = context::MatrixXs::Zero(model.nv, matrix.cols());
        self.jacobianTransposeMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
