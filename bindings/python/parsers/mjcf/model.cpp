//
// Copyright (c) 2024 INRIA
//

#include "pinocchio/parsers/mjcf.hpp"
#include "pinocchio/bindings/python/parsers/mjcf.hpp"
#include "pinocchio/bindings/python/utils/keep-alive.hpp"
#include "pinocchio/bindings/python/utils/path.hpp"

#include <boost/python.hpp>
#include <boost/filesystem/fstream.hpp>

namespace pinocchio
{
  namespace python
  {

    namespace bp = boost::python;

    Model buildModelFromMJCF(const bp::object & filename)
    {
      Model model;
      ::pinocchio::mjcf::buildModel(path(filename), model);
      return model;
    }

    Model & buildModelFromMJCF(const bp::object & filename, Model & model)
    {
      ::pinocchio::mjcf::buildModel(path(filename), model);
      return model;
    }

    Model & buildModelFromMJCF(
      const bp::object & filename,
      const JointModel & root_joint,
      const std::string & root_joint_name,
      Model & model)
    {
      ::pinocchio::mjcf::buildModel(path(filename), root_joint, root_joint_name, model);
      return model;
    }

    Model &
    buildModelFromMJCF(const bp::object & filename, const JointModel & root_joint, Model & model)
    {
      return buildModelFromMJCF(filename, root_joint, "root_joint", model);
    }

    Model buildModelFromMJCF(
      const bp::object & filename,
      const JointModel & root_joint,
      const std::string & root_joint_name = "root_joint")
    {
      Model model;
      return buildModelFromMJCF(filename, root_joint, root_joint_name, model);
    }

    Model buildModelFromMJCFString(const std::string & xml_string)
    {
      Model model;
      ::pinocchio::mjcf::buildModelFromXML(xml_string, model);
      return model;
    }

    Model buildModelFromMJCFString(
      const std::string & xml_string,
      const JointModel & root_joint,
      const std::string & root_joint_name = "root_joint")
    {
      Model model;
      ::pinocchio::mjcf::buildModelFromXML(xml_string, root_joint, root_joint_name, model);
      return model;
    }

    PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel)
    buildBilateralConstraintModelsFromMJCF(Model & model, const bp::object & filename)
    {
      PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel)
      point_anchor_constraint_models;
      ::pinocchio::mjcf::buildConstraintModelsFromXML(
        path(filename), model, point_anchor_constraint_models);
      return point_anchor_constraint_models;
    }

    PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel)
    &buildBilateralConstraintModelsFromMJCF(
      Model & model,
      const bp::object & filename,
      PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel) & point_anchor_constraint_models)
    {
      ::pinocchio::mjcf::buildConstraintModelsFromXML(
        path(filename), model, point_anchor_constraint_models);
      return point_anchor_constraint_models;
    }

    PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel)
    buildFrameAnchorConstraintModelsFromMJCF(Model & model, const bp::object & filename)
    {
      PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel)
      frame_anchor_constraint_models;
      ::pinocchio::mjcf::buildConstraintModelsFromXML(
        path(filename), model, frame_anchor_constraint_models);
      return frame_anchor_constraint_models;
    }

    PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel)
    &buildFrameAnchorConstraintModelsFromMJCF(
      Model & model,
      const bp::object & filename,
      PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel) & frame_anchor_constraint_models)
    {
      ::pinocchio::mjcf::buildConstraintModelsFromXML(
        path(filename), model, frame_anchor_constraint_models);
      return frame_anchor_constraint_models;
    }

    void buildAllConstraintModelsFromMJCF(
      Model & model,
      const bp::object & filename,
      PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel) & point_anchor_constraint_models,
      PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel) & frame_anchor_constraint_models)
    {
      ::pinocchio::mjcf::buildConstraintModelsFromXML(
        path(filename), model, point_anchor_constraint_models, frame_anchor_constraint_models);
    }

    void exposeMJCFModel()
    {
      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &)>(pinocchio::python::buildModelFromMJCF),
        bp::arg("mjcf_filename"),
        "Parse the MJCF file given in input and return a pinocchio Model.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model & (*)(const bp::object &, Model &)>(
          pinocchio::python::buildModelFromMJCF),
        bp::args("mjcf_filename", "model"),
        "Parse the MJCF file given in input and return a pinocchio Model.",
        bp::return_internal_reference<2>());

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &, const JointModel &, const std::string &)>(
          pinocchio::python::buildModelFromMJCF),
        (bp::args("mjcf_filename", "root_joint"), bp::arg("root_joint_name") = "root_joint"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint.");

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model & (*)(const bp::object &, const JointModel &, const std::string &,
                                Model &)>(pinocchio::python::buildModelFromMJCF),
        (bp::args("mjcf_filename", "root_joint"), bp::arg("root_joint_name") = "root_joint",
         bp::arg("model")),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint and the "
        "constraint models.",
        bp::return_internal_reference<4>());

      bp::def(
        "buildModelFromMJCF",
        static_cast<Model (*)(const bp::object &, const JointModel &, const std::string &)>(
          pinocchio::python::buildModelFromMJCF),
        (bp::args("mjcf_filename", "root_joint"), bp::arg("root_joint_name") = "root_joint"),
        "Parse the MJCF file and return a pinocchio Model with the given root Joint and its "
        "specified name as well as a constraint list if some are present in the MJCF file.");

      bp::def(
        "buildModelFromMJCFString",
        static_cast<Model (*)(const std::string &)>(pinocchio::python::buildModelFromMJCFString),
        bp::args("xml_string"),
        "Parse the MJCF string given in input and return a pinocchio Model.");

      bp::def(
        "buildModelFromMJCFString",
        static_cast<Model (*)(const std::string &, const JointModel &, const std::string &)>(
          pinocchio::python::buildModelFromMJCFString),
        (bp::args("mjcf_filename", "root_joint"), bp::arg("root_joint_name") = "root_joint"),
        "Parse the MJCF string and return a pinocchio Model with the given root Joint and its "
        "specified name as well as a constraint list if some are present in the MJCF file.");

      bp::def(
        "buildBilateralConstraintModelsFromMJCF",
        static_cast<PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel) (*)(
          Model &, const bp::object &)>(pinocchio::python::buildBilateralConstraintModelsFromMJCF),
        bp::args("mjcf_filename", "model"),
        "Parse the MJCF file given in input and return a list of pinocchio CosntraintModel.");

      bp::def(
        "buildBilateralConstraintModelsFromMJCF",
        static_cast < PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel)
          & (*)(Model &, const bp::object &,
                PINOCCHIO_ALIGNED_STD_VECTOR(PointAnchorConstraintModel) &)
              > (pinocchio::python::buildBilateralConstraintModelsFromMJCF),
        bp::args("mjcf_filename", "model", "point_anchor_constraint_models"),
        "Parse the MJCF file given in input and return a list of pinocchio CosntraintModel.",
        bp::return_internal_reference<3>());

      bp::def(
        "buildFrameAnchorConstraintModelsFromMJCF",
        static_cast<PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel) (*)(
          Model &, const bp::object &)>(
          pinocchio::python::buildFrameAnchorConstraintModelsFromMJCF),
        bp::args("mjcf_filename", "model"),
        "Parse the MJCF file given in input and return a list of pinocchio CosntraintModel.");

      bp::def(
        "buildFrameAnchorConstraintModelsFromMJCF",
        static_cast < PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel)
          & (*)(Model &, const bp::object &,
                PINOCCHIO_ALIGNED_STD_VECTOR(FrameAnchorConstraintModel) &)
              > (pinocchio::python::buildFrameAnchorConstraintModelsFromMJCF),
        bp::args("mjcf_filename", "model", "frame_anchor_constraint_models"),
        "Parse the MJCF file given in input and return a list of pinocchio CosntraintModel.",
        bp::return_internal_reference<3>());

      bp::def(
        "buildAllConstraintModelsFromMJCF", pinocchio::python::buildAllConstraintModelsFromMJCF,
        bp::args(
          "mjcf_filename", "model", "point_anchor_constraint_models",
          "frame_anchor_constraint_models"),
        "Parse the MJCF file given in input and fill constaint models vectors.");
    }
  } // namespace python
} // namespace pinocchio
