//
// Copyright (c) 2022 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_set_trivial_cones_hpp__
#define __pinocchio_python_algorithm_constraints_set_trivial_cones_hpp__

#include "pinocchio/algorithm/constraints/sets/zero-cone.hpp"
#include "pinocchio/algorithm/constraints/sets/full-space-cone.hpp"
#include "pinocchio/algorithm/constraints/sets/orthant-cone.hpp"

#include "pinocchio/bindings/python/algorithm/constraints/set-base.hpp"
#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename TrivialCone>
    struct TrivialConePythonVisitor
    : public boost::python::def_visitor<TrivialConePythonVisitor<TrivialCone>>
    {

      static void expose(const std::string & class_name, const std::string & doc_string = "")
      {
        bp::class_<TrivialCone>(class_name.c_str(), doc_string.c_str(), bp::no_init)
          .def(SetPythonVisitor<TrivialCone, context::VectorXs>())
          .def(ConeSetPythonVisitor<TrivialCone>())
          .def(CopyableVisitor<TrivialCone>());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_set_trivial_cones_hpp__
