//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_point_constraint_model_base_hpp__
#define __pinocchio_algorithm_constraints_point_constraint_model_base_hpp__

#include <algorithm>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/spatial/skew.hpp"
#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/binary-kinematics-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-common-parameters.hpp"
#include "pinocchio/algorithm/constraints/baumgarte-corrector-parameters.hpp"

namespace pinocchio
{

  template<typename Derived>
  struct PointConstraintModelBase;

  template<typename Derived>
  struct traits<PointConstraintModelBase<Derived>>
  {
    enum
    {
      Size = 3
    };

    // --------------------------------------------------------------
    // Traits characterizing the constraint behaviour in CRTP
    // --------------------------------------------------------------
    static constexpr ConstraintFormulationLevel constraint_formulation_level =
      ConstraintFormulationLevel::VELOCITY_LEVEL;
    static constexpr ConstraintSizeType constraint_size_type = ConstraintSizeType::STATIC;

    static constexpr bool has_baumgarte_corrector =
      true; // Baumgarte make sense and exist directly for the constraint
    static constexpr bool has_compliance_member =
      true; // The constraint itself posses a member m_compliance which can be set by the user
    static constexpr bool has_set = true; // The constraint itself defines the set, otherwise must
                                          // have a mechanism for set-related visitors

    // --------------------------------------------------------------
    // Traits for the algorithmic methods on current state
    // --------------------------------------------------------------
    // Template to generate type
    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::
        Matrix<Scalar, Size, InputMatrixPlain::ColsAtCompileTime, InputMatrixPlain::Options>
          type;
    };

    template<typename InputMatrix>
    struct JacobianTransposeMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Eigen::Dynamic,
        InputMatrixPlain::ColsAtCompileTime,
        InputMatrixPlain::Options>
        type;
    };
  };

  ///
  /// \brief Contact model structure containg all the info describing the rigid contact model
  ///
  template<typename Derived>
  struct PointConstraintModelBase : BinaryKinematicsConstraintModelBase<Derived>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename traits<Derived>::Scalar Scalar;
    static constexpr int Options = traits<Derived>::Options;

    typedef BinaryKinematicsConstraintModelBase<Derived> Base;
    typedef typename Base::BaseCommonParameters BaseCommonParameters;
    typedef ConstraintModelBase<Derived> RootBase;

    template<typename OtherDerived>
    friend struct PointConstraintModelBase;

    static const ConstraintFormulationLevel constraint_formulation_level =
      traits<PointConstraintModelBase>::constraint_formulation_level;
    typedef typename traits<Derived>::ConstraintData ConstraintData;
    typedef typename traits<Derived>::ComplianceVectorType ComplianceVectorType;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef Eigen::Matrix<Scalar, 3, 6, Options> Matrix36;
    typedef Eigen::Matrix<Scalar, 6, 6, Options> Matrix6;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;
    typedef Vector3 VectorConstraintSize;

    using Base::getA1;
    using Base::getA2;
    using RootBase::jacobianMatrixProduct;
    using RootBase::jacobianTransposeMatrixProduct;
    using RootBase::residualSize;

    // -------------------------------
    // METHODS SPECIFIC TO CLASS
    // -------------------------------

    // CRTP related ------------------

    /// \brief Cast to Base
    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    /// \brief Const cast to Base
    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    /// \brief Cast to BaseCommonParameters
    BaseCommonParameters & base_common_parameters()
    {
      return static_cast<BaseCommonParameters &>(*this);
    }

    /// \brief Const cast to BaseCommonParameters
    const BaseCommonParameters & base_common_parameters() const
    {
      return static_cast<const BaseCommonParameters &>(*this);
    }

    // Constructors ------------------

    ///
    /// \brief Default constructor
    ///
    PointConstraintModelBase()
    : Base()
    {
    }

    ///
    /// \brief Contructor with from a given type, joint indexes and placements.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] model Model associated to the constraint.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// \param[in] joint2_placement Placement of the constraint w.r.t the frame of joint2.
    /// \param[in] reference_frame Reference frame in which the constraints quantities are
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement)
    : Base(model, joint1_id, joint1_placement, joint2_id, joint2_placement)
    {
    }

    ///
    /// \brief Contructor with from a given type and .
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    ///
    /// \remarks The second joint id (joint2_id) is set to be 0 (corresponding to the index of the
    /// universe).
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model, const JointIndex joint1_id)
    : Base(model, joint1_id)
    {
    }

    ///
    /// \brief Contructor with from a given type, joint1_id and placement.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint1_placement Placement of the constraint w.r.t the frame of joint1.
    /// expressed.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement)
    : Base(model, joint1_id, joint1_placement)
    {
    }

    ///
    /// \brief Contructor with from a given type and the joint ids.
    ///
    /// \param[in] type Type of the contact.
    /// \param[in] joint1_id Index of the joint 1 in the model tree.
    /// \param[in] joint2_id Index of the joint 2 in the model tree.
    ///
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointConstraintModelBase(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id)
    : Base(model, joint1_id, joint2_id)
    {
    }

    // Operators ---------------------

    ///
    /// \brief Comparison operator
    ///
    /// \param[in] other Other PointConstraintModelBase to compare with.
    ///
    /// \returns true if the two *this is equal to other (type, joint1_id and placement attributs
    /// must be the same).
    ///
    bool operator==(const PointConstraintModelBase & other) const
    {
      return Base::operator==(other);
    }

    ///
    /// \brief Opposite of the comparison operator.
    ///
    /// \param[in] other Other PointConstraintModelBase to compare with.
    ///
    /// \returns false if the two *this is not equal to other (at least type, joint1_id or placement
    /// attributs is different).
    ///
    bool operator!=(const PointConstraintModelBase & other) const
    {
      return Base::operator!=(other);
    }

    /// \returns Cast to NewScalar
    template<typename NewScalar, typename OtherDerived>
    void cast(PointConstraintModelBase<OtherDerived> & res) const
    {
      Base::template cast<NewScalar>(res);
    }

    // Methods for rigid body --------

    ///
    /// \brief This function computes the spatial inertia associated with the constraint.
    /// This function is useful to express the constraint inertia associated with the constraint for
    /// AL-based approaches.
    ///
    template<typename Vector3Like>
    Matrix6 computeConstraintSpatialInertia(
      const SE3Tpl<Scalar, Options> & placement,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia) const
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector3Like, Vector3);
      Matrix6 res;

      const auto & R = placement.rotation();
      const auto & t = placement.translation();

      typedef Eigen::Matrix<Scalar, 3, 3, Options> Matrix3;
      const Matrix3 R_Sigma = R * diagonal_constraint_inertia.asDiagonal();
      const Matrix3 t_skew = skew(t);

      auto block_LL = res.template block<3, 3>(SE3::LINEAR, SE3::LINEAR);
      auto block_LA = res.template block<3, 3>(SE3::LINEAR, SE3::ANGULAR);
      auto block_AL = res.template block<3, 3>(SE3::ANGULAR, SE3::LINEAR);
      auto block_AA = res.template block<3, 3>(SE3::ANGULAR, SE3::ANGULAR);

      block_LL.noalias() = R_Sigma * R.transpose();
      block_LA.noalias() = -block_LL * t_skew;
      block_AL.noalias() = block_LA.transpose();
      block_AA.noalias() = t_skew * block_LA;

      return res;
    }

    template<
      typename Matrix6LikeOut1,
      typename Matrix6LikeOut2,
      typename Matrix6LikeOut3,
      ReferenceFrame rf>
    void computeConstraintInertias(
      const ConstraintData & cdata,
      const Scalar & constraint_inertia_value,
      const Eigen::MatrixBase<Matrix6LikeOut1> & I11,
      const Eigen::MatrixBase<Matrix6LikeOut2> & I12,
      const Eigen::MatrixBase<Matrix6LikeOut3> & I22,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      computeConstraintInertias(
        cdata, Vector3::Constant(constraint_inertia_value), I11.const_cast_derived(),
        I12.const_cast_derived(), I22.const_cast_derived(), reference_frame);
    }

    template<
      typename Vector3Like,
      typename Matrix6LikeOut1,
      typename Matrix6LikeOut2,
      typename Matrix6LikeOut3,
      ReferenceFrame rf>
    void computeConstraintInertias(
      const ConstraintData & cdata,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia,
      const Eigen::MatrixBase<Matrix6LikeOut1> & I11,
      const Eigen::MatrixBase<Matrix6LikeOut2> & I12,
      const Eigen::MatrixBase<Matrix6LikeOut3> & I22,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector3Like, Vector3);
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Matrix6LikeOut1, Matrix6);
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Matrix6LikeOut2, Matrix6);
      EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Matrix6LikeOut3, Matrix6);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);
      //      assert((check_expression_if_real<Scalar,
      //      true>(diagonal_constraint_inertia.isZero(Scalar(0)))));

      const auto & A1 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A1_world : cdata.A1_local;
      const auto & A2 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A2_world : cdata.A2_local;

      Matrix36 diagonal_constraint_inertia_time_A;

      if (this->joint1_id > 0)
      {
        diagonal_constraint_inertia_time_A.noalias() =
          diagonal_constraint_inertia.asDiagonal() * A1;
        I11.const_cast_derived().noalias() = A1.transpose() * diagonal_constraint_inertia_time_A;
      }
      else
        I11.const_cast_derived().setZero();

      if (this->joint2_id > 0)
      {
        diagonal_constraint_inertia_time_A.noalias() =
          diagonal_constraint_inertia.asDiagonal() * A2;
        I22.const_cast_derived().noalias() = A2.transpose() * diagonal_constraint_inertia_time_A;
      }
      else
        I22.const_cast_derived().setZero();

      // Compute the cross coupling term
      if (this->joint1_id > 0 && this->joint2_id > 0)
      {
        I12.const_cast_derived().noalias() = A1.transpose() * diagonal_constraint_inertia_time_A;
      }
      else
        I12.const_cast_derived().setZero();
    }

    // -------------------------------
    // IMPLEMENTATIONS OF BASE METHODS
    // -------------------------------

    // Methods for algorithms --------

    /// \copydoc RootBase::calc
    template<template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);

      if (this->joint1_id > 0)
        cdata.oMc1 = data.oMi[this->joint1_id] * this->joint1_placement;
      else
        cdata.oMc1 = this->joint1_placement;

      if (this->joint2_id > 0)
        cdata.oMc2 = data.oMi[this->joint2_id] * this->joint2_placement;
      else
        cdata.oMc2 = this->joint2_placement;

      // Compute relative placement
      cdata.c1Mc2 = cdata.oMc1.actInv(cdata.oMc2);
      const auto & _1R2_ = cdata.c1Mc2.rotation();

      // Compute errors
      auto & position_error = cdata.constraint_position_error;
      position_error.noalias() = cdata.c1Mc2.translation();
      position_error.noalias() -= this->desired_constraint_offset;
      //      cdata.constraint_position_error = cdata.oMc1.inverse().translation();

      const auto vf1 = this->joint1_placement.actInv(data.v[this->joint1_id]);
      const auto vf2 = this->joint2_placement.actInv(data.v[this->joint2_id]);

      auto & velocity_error = cdata.constraint_velocity_error;
      const Vector3 velocity_error_component1 = _1R2_ * vf2.linear() - vf1.linear();
      velocity_error.noalias() = velocity_error_component1 - vf1.angular().cross(position_error);
      velocity_error.noalias() -= this->desired_constraint_velocity;

      const auto af1 = this->joint1_placement.actInv(data.a[this->joint1_id]);
      const auto af2 = this->joint2_placement.actInv(data.a[this->joint2_id]);
      auto & acceleration_error = cdata.constraint_acceleration_error;
      acceleration_error.noalias() = _1R2_ * (af2.linear() + vf2.angular().cross(vf2.linear()))
                                     - (af1.linear() + vf1.angular().cross(vf1.linear()));
      acceleration_error.noalias() -= af1.angular().cross(position_error);
      acceleration_error.noalias() += vf1.angular().cross(vf1.angular().cross(position_error));
      acceleration_error.noalias() -= Scalar(2) * vf1.angular().cross(velocity_error_component1);
      acceleration_error.noalias() -= this->desired_constraint_acceleration;

      cdata.A1_world = this->getA1(cdata, WorldFrameTag());
      cdata.A2_world = this->getA2(cdata, WorldFrameTag());
      cdata.A_world = cdata.A1_world + cdata.A2_world;

      cdata.A1_local = this->getA1(cdata, LocalFrameTag());
      cdata.A2_local = this->getA2(cdata, LocalFrameTag());
      cdata.A_local = cdata.A1_local + cdata.A2_local;
    }

    /// \copydoc RootBase::jacobian
    template<template<typename, int> class JointCollectionTpl, typename JacobianMatrix>
    void jacobianImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

      const SE3 & oMc1 = cdata.oMc1;
      const SE3 & oMc2 = cdata.oMc2;
      const SE3 & c1Mc2 = cdata.c1Mc2;
      const auto & position_error = cdata.constraint_position_error;

      for (Eigen::Index jj = 0; jj < model.nv; ++jj)
      {
        if (this->colwise_joint1_sparsity[jj] || this->colwise_joint2_sparsity[jj])
        {
          typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
          const ConstColXpr Jcol = data.J.col(jj);
          const MotionRef<const ConstColXpr> Jcol_motion(Jcol);

          jacobian_matrix.col(jj).setZero();
          if (this->colwise_joint1_sparsity[jj])
          {
            const Motion Jcol_local(oMc1.actInv(Jcol_motion)); // TODO: simplify computations
            jacobian_matrix.col(jj).noalias() -= Jcol_local.linear();
            jacobian_matrix.col(jj).noalias() += -Jcol_local.angular().cross(position_error);
          }

          if (this->colwise_joint2_sparsity[jj])
          {
            const Motion Jcol_local(oMc2.actInv(Jcol_motion)); // TODO: simplify computations
            jacobian_matrix.col(jj) += c1Mc2.rotation() * Jcol_local.linear();
          }
        }
      }
    }

    /// \copydoc RootBase::jacobianMatrixProduct
    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Derived>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(3, mat.cols());
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    /// \copydoc RootBase::jacobianMatrixProduct
    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef typename Data::Vector3 Vector3;
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), residualSize(cdata));
      PINOCCHIO_UNUSED_VARIABLE(aot);

      if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      //      const Eigen::Index constraint_size = size();
      //
      //      const Eigen::Index
      //      complexity_strategy_1 = 6 * res.cols() * 36 + constraint_size * 36 * res.cols(),
      //      complexity_strategy_2 = 36 * constraint_size * 6 + constraint_size * 36 * res.cols();

      const auto & A1 = cdata.A1_world;
      const auto & A2 = cdata.A2_world;

      const auto & A = cdata.A_world;
      for (Eigen::Index jj = 0; jj < model.nv; ++jj)
      {
        if (!(this->colwise_joint1_sparsity[jj] || this->colwise_joint2_sparsity[jj]))
          continue;
        Vector3 AxSi;

        typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);

        if (this->colwise_joint1_sparsity[jj] && this->colwise_joint2_sparsity[jj])
        {
          AxSi.noalias() = A * Jcol;
        }
        else if (this->colwise_joint1_sparsity[jj])
          AxSi.noalias() = A1 * Jcol;
        else
          AxSi.noalias() = A2 * Jcol;

        if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.noalias() -= AxSi * mat.row(jj);
        else // AddTo, SetTo
          res.noalias() += AxSi * mat.row(jj);
      }
    }

    /// \copydoc RootBase::jacobianTransposeMatrixProduct
    template<typename InputMatrix, template<typename, int> class JointCollectionTpl>
    typename traits<Derived>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Derived>::template JacobianTransposeMatrixProductReturnType<
        InputMatrix>::type ReturnType;
      ReturnType res(model.nv, mat.cols());
      jacobianTransposeMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    /// \copydoc RootBase::jacobianTransposeMatrixProduct
    template<
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op = SETTO>
    void jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot = SetTo()) const
    {
      typedef DataTpl<Scalar, Options, JointCollectionTpl> Data;
      typedef typename Data::Vector3 Vector3;
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), residualSize(cdata));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);
      PINOCCHIO_UNUSED_VARIABLE(aot);

      if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      const auto & A1 = cdata.A1_world;
      const auto & A2 = cdata.A2_world;

      const auto & A = cdata.A_world;
      for (Eigen::Index jj = 0; jj < model.nv; ++jj)
      {
        if (!(this->colwise_joint1_sparsity[jj] || this->colwise_joint2_sparsity[jj]))
          continue;
        Vector3 AxSi;

        typedef typename Data::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);

        if (this->colwise_joint1_sparsity[jj] && this->colwise_joint2_sparsity[jj])
        {
          AxSi.noalias() = A * Jcol;
        }
        else if (this->colwise_joint1_sparsity[jj])
          AxSi.noalias() = A1 * Jcol;
        else
          AxSi.noalias() = A2 * Jcol;

        if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.row(jj).noalias() -= AxSi.transpose() * mat;
        else
          res.row(jj).noalias() += AxSi.transpose() * mat;
      }
    }

    /// \copydoc RootBase::mapConstraintForceToJointForces
    template<
      template<typename, int> class JointCollectionTpl,
      typename ForceLike,
      typename ForceAllocator,
      ReferenceFrame rf>
    void mapConstraintForceToJointForcesImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, Options>, ForceAllocator> & joint_forces,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), residualSize(cdata));
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      // Todo: optimize code
      const auto & A1 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A1_world : cdata.A1_local;
      const auto & A2 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A2_world : cdata.A2_local;

      if (this->joint1_id > 0)
        joint_forces[this->joint1_id].toVector().noalias() +=
          A1.transpose() * constraint_forces.template head<3>();
      if (this->joint2_id > 0)
        joint_forces[this->joint2_id].toVector().noalias() +=
          A2.transpose() * constraint_forces.template head<3>();
    }

    /// \copydoc RootBase::mapJointMotionsToConstraintMotion
    template<
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointMotionsToConstraintMotionImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      const DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, Options>, MotionAllocator> & joint_accelerations,
      const Eigen::MatrixBase<VectorLike> & constraint_motion,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_accelerations.size(), size_t(model.njoints));
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motion.rows(), residualSize(cdata));
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      const auto & A1 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A1_world : cdata.A1_local;
      const auto & A2 =
        std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value ? cdata.A2_world : cdata.A2_local;

      if (this->joint1_id > 0 && this->joint2_id > 0)
        constraint_motion.const_cast_derived().template head<3>().noalias() =
          A1 * joint_accelerations[this->joint1_id].toVector()
          + A2 * joint_accelerations[this->joint2_id].toVector();
      else if (this->joint1_id > 0)
        constraint_motion.const_cast_derived().template head<3>().noalias() =
          A1 * joint_accelerations[this->joint1_id].toVector();
      else if (this->joint2_id > 0)
        constraint_motion.const_cast_derived().template head<3>().noalias() =
          A2 * joint_accelerations[this->joint2_id].toVector();
      else
        constraint_motion.const_cast_derived().setZero();
    }

    /// \copydoc RootBase::appendCouplingConstraintInertias
    template<
      template<typename, int> class JointCollectionTpl,
      typename Vector3Like,
      ReferenceFrame rf>
    void appendCouplingConstraintInertiasImpl(
      const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
      DataTpl<Scalar, Options, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<Vector3Like> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);

      Matrix6 I11, I12, I22;
      computeConstraintInertias(cdata, diagonal_constraint_inertia, I11, I12, I22, reference_frame);
      assert(
        (std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value
         || std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
        && "must never happened");

      Matrix6 & Y1 = std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value
                       ? data.oYaba_augmented[this->joint1_id]
                       : data.oYaba_augmented[this->joint1_id];

      if (this->joint1_id > 0)
        Y1 += I11;

      Matrix6 & Y2 = std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value
                       ? data.oYaba_augmented[this->joint2_id]
                       : data.oYaba_augmented[this->joint2_id];

      if (this->joint2_id > 0)
        Y2 += I22;

      if (this->joint1_id > 0 && this->joint2_id > 0)
      {
        assert(
          data.joint_cross_coupling.exists({this->joint1_id, this->joint2_id})
          || data.joint_cross_coupling.exists({this->joint2_id, this->joint1_id}));
        if (data.joint_cross_coupling.exists({this->joint1_id, this->joint2_id}))
        {
          data.joint_cross_coupling.get({this->joint1_id, this->joint2_id}) += I12;
        }
        else
        {
          data.joint_cross_coupling.get({this->joint2_id, this->joint1_id}) += I12.transpose();
        }
      }
    }

    /// \copydoc Base::getA1
    template<ReferenceFrame rf>
    Matrix36 getA1Impl(const ConstraintData & cdata, ReferenceFrameTag<rf>) const
    {
      Matrix36 res;

      if constexpr (std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = oM1.rotation().transpose() * v_tmp;

        const SE3 & oM1 = cdata.oMc1;
        Vector3 v_tmp;
        res.template leftCols<3>() = -oM1.rotation().transpose();
        INTERNAL_LOOP(0, -oM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, -oM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, -oM1.translation(), res.template rightCols<3>());

        for (int i = 0; i < 3; ++i)
        {
          res.template rightCols<3>().col(i) +=
            cdata.constraint_position_error.cross(oM1.rotation().transpose().col(i));
        }

#undef INTERNAL_LOOP
      }
      else if constexpr (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = iM1.rotation().transpose() * v_tmp;

        const SE3 & iM1 = this->joint1_placement;
        Vector3 v_tmp;
        res.template leftCols<3>() = -iM1.rotation().transpose();
        INTERNAL_LOOP(0, -iM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, -iM1.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, -iM1.translation(), res.template rightCols<3>());

        for (int i = 0; i < 3; ++i)
        {
          res.template rightCols<3>().col(i) +=
            cdata.constraint_position_error.cross(iM1.rotation().transpose().col(i));
        }
      }
      else
      {
        assert(false && "Should never happened");
      }

#undef INTERNAL_LOOP

      return res;
    }

    /// \copydoc Base::getA2
    template<ReferenceFrame rf>
    Matrix36 getA2Impl(const ConstraintData & cdata, ReferenceFrameTag<rf>) const
    {
      Matrix36 res;
      typedef typename SE3::Vector3 Vector3;

      if constexpr (std::is_same<ReferenceFrameTag<rf>, WorldFrameTag>::value)
      {
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = oM1.rotation().transpose() * v_tmp;

        const SE3 & oM1 = cdata.oMc1;
        const SE3 & oM2 = cdata.oMc2;
        res.template leftCols<3>() = oM1.rotation().transpose();
        Vector3 v_tmp;
        INTERNAL_LOOP(0, oM2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, oM2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, oM2.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }
      else if constexpr (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
      {
        const SE3 & j2Mc2 = this->joint2_placement;
        const SE3 & c1Mc2 = cdata.c1Mc2;
        const typename SE3::Matrix3 c1Rj2 = c1Mc2.rotation() * j2Mc2.rotation().transpose();
        res.template leftCols<3>() = c1Rj2;
        Vector3 v_tmp;
#define INTERNAL_LOOP(axis_id, v3_in, res)                                                         \
  CartesianAxis<axis_id>::cross(v3_in, v_tmp);                                                     \
  res.col(axis_id).noalias() = c1Rj2 * v_tmp;

        INTERNAL_LOOP(0, j2Mc2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(1, j2Mc2.translation(), res.template rightCols<3>());
        INTERNAL_LOOP(2, j2Mc2.translation(), res.template rightCols<3>());

#undef INTERNAL_LOOP
      }
      else
      {
        assert(false && "Should never happened");
      }

      return res;
    }
  }; // PointConstraintModelBase<Derived>

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraints_point_constraint_model_base_hpp__
