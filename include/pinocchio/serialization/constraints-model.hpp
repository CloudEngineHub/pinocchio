//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_constraints_model_hpp__
#define __pinocchio_serialization_constraints_model_hpp__

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/serialization/eigen.hpp"
#include "pinocchio/serialization/se3.hpp"
#include "pinocchio/serialization/eigen-storage.hpp"
#include "pinocchio/serialization/constraints-set.hpp"
#include "pinocchio/serialization/boost-blank.hpp"

#include <boost/serialization/variant.hpp>

namespace boost
{
  namespace serialization
  {

    template<typename Archive, typename Scalar>
    void serialize(
      Archive & ar,
      ::pinocchio::BaumgarteCorrectorParametersTpl<Scalar> & baumgarte_parameters,
      const unsigned int /*version*/)
    {
      ar & make_nvp("Kp", baumgarte_parameters.Kp);
      ar & make_nvp("Kd", baumgarte_parameters.Kd);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar, ::pinocchio::ConstraintModelBase<Derived> & cmodel, const unsigned int version)
    {
      PINOCCHIO_UNUSED_VARIABLE(version);
      ar & make_nvp("name", cmodel.name);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::KinematicsConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::KinematicsConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::JointWiseConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointWiseConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    namespace internal
    {
      template<typename Derived>
      struct ConstraintModelCommonParametersAccessor
      : public ::pinocchio::ConstraintModelCommonParameters<Derived>
      {
        typedef ::pinocchio::ConstraintModelCommonParameters<Derived> Base;
        using Base::m_baumgarte_parameters;
        using Base::m_compliance;
      };
    } // namespace internal

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::ConstraintModelCommonParameters<Derived> & cmodel,
      const unsigned int version)
    {
      PINOCCHIO_UNUSED_VARIABLE(version);
      typedef internal::ConstraintModelCommonParametersAccessor<Derived> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("m_compliance", cmodel_.m_compliance);
      ar & make_nvp("m_baumgarte_parameters", cmodel_.m_baumgarte_parameters);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct JointLimitConstraintModelAccessor
      : public ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options> Base;

        using Base::activable_idx_in_selected;
        using Base::activable_idx_qs;
        using Base::activable_idx_qs_reduce;
        using Base::activable_position_limit;
        using Base::activable_position_margin;
        using Base::lower_max_residual_size;
        using Base::max_of_nvs;
        using Base::nq_reduce;
        using Base::selected_joint_idx_vs;
        using Base::selected_joint_nqs;
        using Base::selected_joint_nvs;
        using Base::selected_joints;
        using Base::selected_row_indexes;
        using Base::selected_row_sparsity_pattern;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointLimitConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      typedef internal::JointLimitConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("selected_joints", cmodel_.selected_joints);
      ar & make_nvp("selected_row_sparsity_pattern", cmodel_.selected_row_sparsity_pattern);
      ar & make_nvp("selected_row_indexes", cmodel_.selected_row_indexes);
      ar & make_nvp("selected_joint_nqs", cmodel_.selected_joint_nqs);
      ar & make_nvp("selected_joint_nvs", cmodel_.selected_joint_nvs);
      ar & make_nvp("selected_joint_idx_vs", cmodel_.selected_joint_idx_vs);
      ar & make_nvp("nq_reduce", cmodel_.nq_reduce);
      ar & make_nvp("max_of_nvs", cmodel_.max_of_nvs);
      ar & make_nvp("activable_idx_in_selected", cmodel_.activable_idx_in_selected);
      ar & make_nvp("activable_idx_qs", cmodel_.activable_idx_qs);
      ar & make_nvp("activable_idx_qs_reduce", cmodel_.activable_idx_qs_reduce);
      ar & make_nvp("activable_position_limit", cmodel_.activable_position_limit);
      ar & make_nvp("activable_position_margin", cmodel_.activable_position_margin);
      ar & make_nvp("lower_max_residual_size", cmodel_.lower_max_residual_size);
    }

    namespace internal
    {
      template<typename Scalar, int Options>
      struct JointFrictionConstraintModelAccessor
      : public ::pinocchio::JointFrictionConstraintModelTpl<Scalar, Options>
      {
        typedef ::pinocchio::JointFrictionConstraintModelTpl<Scalar, Options> Base;
        using Base::active_dofs;
        using Base::active_joints;
        using Base::m_friction_lower_limit;
        using Base::m_friction_upper_limit;
        using Base::row_active_indexes;
        using Base::row_sparsity_pattern;
      };
    } // namespace internal

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::JointFrictionConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::JointFrictionConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      typedef internal::JointFrictionConstraintModelAccessor<Scalar, Options> Accessor;
      auto & cmodel_ = reinterpret_cast<Accessor &>(cmodel);
      ar & make_nvp("active_joints", cmodel_.active_joints);
      ar & make_nvp("active_dofs", cmodel_.active_dofs);
      ar & make_nvp("row_sparsity_pattern", cmodel_.row_sparsity_pattern);
      ar & make_nvp("row_active_indexes", cmodel_.row_active_indexes);
      ar & make_nvp("friction_lower_limit", cmodel_.m_friction_lower_limit);
      ar & make_nvp("friction_upper_limit", cmodel_.m_friction_upper_limit);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::BinaryKinematicsConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::BinaryKinematicsConstraintModelBase<Derived> Self;
      typedef typename Self::Base KinematicsBase;
      ar & make_nvp("base", boost::serialization::base_object<KinematicsBase>(cmodel));
      typedef typename Self::BaseCommonParameters BaseCommonParameters;
      ar & make_nvp(
        "base_common_parameters", boost::serialization::base_object<BaseCommonParameters>(cmodel));

      ar & make_nvp("joint1_id", cmodel.joint1_id);
      ar & make_nvp("joint2_id", cmodel.joint2_id);
      ar & make_nvp("joint1_placement", cmodel.joint1_placement);
      ar & make_nvp("joint2_placement", cmodel.joint2_placement);
      ar & make_nvp("desired_constraint_offset", cmodel.desired_constraint_offset);
      ar & make_nvp("desired_constraint_velocity", cmodel.desired_constraint_velocity);
      ar & make_nvp("desired_constraint_acceleration", cmodel.desired_constraint_acceleration);
      ar & make_nvp("colwise_joint1_sparsity", cmodel.colwise_joint1_sparsity);
      ar & make_nvp("colwise_joint2_sparsity", cmodel.colwise_joint2_sparsity);
      ar & make_nvp("joint1_span_indexes", cmodel.joint1_span_indexes);
      ar & make_nvp("joint2_span_indexes", cmodel.joint2_span_indexes);
      ar & make_nvp("colwise_sparsity", cmodel.colwise_sparsity);
      ar & make_nvp("colwise_span_indexes", cmodel.colwise_span_indexes);
      ar & make_nvp("nv", cmodel.nv);
      ar & make_nvp("depth_joint1", cmodel.depth_joint1);
      ar & make_nvp("depth_joint2", cmodel.depth_joint2);
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::PointConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::PointConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::PointAnchorConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::PointAnchorConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::PointContactModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::PointContactModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    template<typename Archive, typename Derived>
    void serialize(
      Archive & ar,
      ::pinocchio::FrameConstraintModelBase<Derived> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrameConstraintModelBase<Derived> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    template<typename Archive, typename Scalar, int Options>
    void serialize(
      Archive & ar,
      ::pinocchio::FrameAnchorConstraintModelTpl<Scalar, Options> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::FrameAnchorConstraintModelTpl<Scalar, Options> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));
    }

    template<
      typename Archive,
      typename Scalar,
      int Options,
      template<typename, int> class ConstraintCollectionTpl>
    void serialize(
      Archive & ar,
      pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> & cmodel,
      const unsigned int /*version*/)
    {
      typedef ::pinocchio::ConstraintModelTpl<Scalar, Options, ConstraintCollectionTpl> Self;
      typedef typename Self::Base Base;
      ar & make_nvp("base", boost::serialization::base_object<Base>(cmodel));

      typedef typename ConstraintCollectionTpl<Scalar, Options>::ConstraintModelVariant
        ConstraintModelVariant;
      ar & make_nvp(
        "base_variant", boost::serialization::base_object<ConstraintModelVariant>(cmodel));
    }

  } // namespace serialization
} // namespace boost

#endif // __pinocchio_serialization_constraints_model_hpp__
