//
// Copyright (c) 2019-2020 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/serialization.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/serialization.hpp"
#endif // PINOCCHIO_LSP

namespace boost
{
  namespace serialization
  {
    template<
      class Archive,
      typename Scalar,
      int Options,
      template<typename, int> class JointCollectionTpl>
    void serialize(
      Archive & ar,
      pinocchio::ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const unsigned int /*version*/)
    {
      ar & make_nvp("nq", model.nq);
      ar & make_nvp("nqs", model.nqs);
      ar & make_nvp("idx_qs", model.idx_qs);
      ar & make_nvp("nv", model.nv);
      ar & make_nvp("nvs", model.nvs);
      ar & make_nvp("idx_vs", model.idx_vs);
      ar & make_nvp("nvExtended", model.nvExtended);
      ar & make_nvp("nvExtendeds", model.nvExtendeds);
      ar & make_nvp("idx_vExtendeds", model.idx_vExtendeds);
      ar & make_nvp("njoints", model.njoints);
      ar & make_nvp("nbodies", model.nbodies);
      ar & make_nvp("nframes", model.nframes);
      ar & make_nvp("parents", model.parents);
      ar & make_nvp("children", model.children);
      ar & make_nvp("names", model.names);
      ar & make_nvp("supports", model.supports);
      ar & make_nvp("mimic_joint_supports", model.mimic_joint_supports);
      ar & make_nvp("subtrees", model.subtrees);
      ar & make_nvp("mimicking_joints", model.mimicking_joints);
      ar & make_nvp("mimicked_joints", model.mimicked_joints);
      ar & make_nvp("gravity", model.gravity);
      ar & make_nvp("name", model.name);

      ar & make_nvp("referenceConfigurations", model.referenceConfigurations);
      ar & make_nvp("armature", model.armature);
      ar & make_nvp("rotorInertia", model.rotorInertia);
      ar & make_nvp("rotorGearRatio", model.rotorGearRatio);
      ar & make_nvp("lowerDryFrictionLimit", model.lowerDryFrictionLimit);
      ar & make_nvp("upperDryFrictionLimit", model.upperDryFrictionLimit);
      ar & make_nvp("damping", model.damping);
      ar & make_nvp("lowerEffortLimit", model.lowerEffortLimit);
      ar & make_nvp("upperEffortLimit", model.upperEffortLimit);
      ar & make_nvp("lowerVelocityLimit", model.lowerVelocityLimit);
      ar & make_nvp("upperVelocityLimit", model.upperVelocityLimit);
      ar & make_nvp("lowerPositionLimit", model.lowerPositionLimit);
      ar & make_nvp("upperPositionLimit", model.upperPositionLimit);
      ar & make_nvp("positionLimitMargin", model.positionLimitMargin);

      ar & make_nvp("inertias", model.inertias);
      ar & make_nvp("jointPlacements", model.jointPlacements);

      ar & make_nvp("joints", model.joints);
      ar & make_nvp("frames", model.frames);
    }

  } // namespace serialization
} // namespace boost
