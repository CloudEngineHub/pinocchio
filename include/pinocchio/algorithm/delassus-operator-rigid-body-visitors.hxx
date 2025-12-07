//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__

#include "pinocchio/utils/eigen.hpp"

#define PROMOTE_STATIC_EVAL(expression) promote_static_eval<0>(expression)
// #define PROMOTE_STATIC_EVAL(expression) expression
#define DO_NOT_PROMOTE_STATIC_EVAL(expression) expression

namespace pinocchio
{

  template<typename DelassusOperator, bool apply_on_the_right = true, bool solve_in_place = true>
  struct DelassusOperatorRigidBodySystemsComputeBackwardPass
  : public fusion::JointUnaryVisitorBase<DelassusOperatorRigidBodySystemsComputeBackwardPass<
      DelassusOperator,
      apply_on_the_right,
      solve_in_place>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename Model::Scalar Scalar;

    typedef boost::fusion::vector<const Model &, Data &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Matrix6 Matrix6;
      typedef typename JointModel::JointDataDerived JointData;
      typedef std::pair<JointIndex, JointIndex> JointPair;

      const JointIndex joint_i = jmodel.id();
      const JointIndex parent = model.parents[joint_i];

      // ApplyOnTheRight
      if (apply_on_the_right)
      {
        auto & Ia = data.Yaba[joint_i];
        jmodel.calc_aba(
          jdata.derived(), jmodel.jointVelocitySelector(model.armature), Ia, parent > 0);
        if (parent > 0)
        {
          data.Yaba[parent] += impl::internal::SE3actOn<Scalar>::run(data.liMi[joint_i], Ia);
        }
      }

      // SolveInPlace
      if (solve_in_place)
      {
        JointData & _jdata_augmented = boost::get<JointData>(data.joints_augmented[joint_i]);
        JointDataBase<JointData> & jdata_augmented =
          static_cast<JointDataBase<JointData> &>(_jdata_augmented);

        auto Jcols = jmodel.jointCols(data.J);
        auto & Ia_augmented = data.oYaba_augmented[joint_i];

        DO_NOT_PROMOTE_STATIC_EVAL(jdata_augmented.U().noalias()) = Ia_augmented * Jcols;
        DO_NOT_PROMOTE_STATIC_EVAL(jdata_augmented.StU().noalias()) =
          Jcols.transpose() * jdata_augmented.U();

        // Account for the rotor inertia contribution
        jdata_augmented.StU() += data.joint_apparent_inertia[joint_i];

        ::pinocchio::matrix_inversion(jdata_augmented.StU(), jdata_augmented.Dinv());

        DO_NOT_PROMOTE_STATIC_EVAL(jdata_augmented.UDinv().noalias()) =
          jdata_augmented.U() * jdata_augmented.Dinv();

        if (parent > 0)
        {
          DO_NOT_PROMOTE_STATIC_EVAL(Ia_augmented.noalias()) -=
            jdata_augmented.UDinv() * jdata_augmented.U().transpose();
          data.oYaba_augmented[parent] += Ia_augmented;
        }

        // End of the classic ABA backward pass - beginning of cross-coupling handling
        const auto & neighbours = data.joint_neighbours;
        auto & joint_cross_coupling = data.joint_cross_coupling;
        auto & projected_joint_cross_coupling = data.projected_joint_cross_coupling;
        const auto & joint_neighbours = neighbours[joint_i];

        if (joint_neighbours.size() == 0)
          return; // We can return from this point as this joint has no neighbours

        using Matrix6xNV = typename std::remove_reference<typename JointData::UDTypeRef>::type;
        typedef Eigen::Map<Matrix6xNV, EIGEN_DEFAULT_ALIGN_BYTES> MapMatrix6xNV;
        MapMatrix6xNV mat1_tmp = MapMatrix6xNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, 6, jmodel.nv()));
        MapMatrix6xNV mat2_tmp = MapMatrix6xNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, 6, jmodel.nv()));

        auto & JDinv = mat1_tmp;
        DO_NOT_PROMOTE_STATIC_EVAL(JDinv.noalias()) = Jcols * jdata_augmented.Dinv();

        // oL == data.oL[joint_i]
        Matrix6 oL = -JDinv * jdata_augmented.U().transpose();
        oL += Matrix6::Identity();

        for (size_t j = 0; j < joint_neighbours.size(); j++)
        {
          const JointIndex joint_j = joint_neighbours[j];

          assert(joint_cross_coupling.exists(JointPair(joint_j, joint_i)));
          const auto & crosscoupling_ji = joint_cross_coupling.get(JointPair(joint_j, joint_i));
          // assert(crosscoupling_ji.isApprox(crosscoupling_ji2));

          assert(projected_joint_cross_coupling.exists(JointPair(joint_j, joint_i)));
          auto & crosscoupling_ji_Jcols =
            projected_joint_cross_coupling[JointPair(joint_j, joint_i)];

          DO_NOT_PROMOTE_STATIC_EVAL(crosscoupling_ji_Jcols.noalias()) =
            crosscoupling_ji * Jcols; // Warning: UDinv() is actually edge_ij * J

          static_assert(
            !PINOCCHIO_DECLTYPE(crosscoupling_ji_Jcols)::IsRowMajor
              && !PINOCCHIO_DECLTYPE(crosscoupling_ji)::IsRowMajor
              && !PINOCCHIO_DECLTYPE(Jcols)::IsRowMajor,
            "All the elements should be of same.");

          auto & crosscoupling_ji_Jcols_Dinv = mat2_tmp;
          DO_NOT_PROMOTE_STATIC_EVAL(crosscoupling_ji_Jcols_Dinv.noalias()) =
            crosscoupling_ji_Jcols * jdata_augmented.Dinv();

          DO_NOT_PROMOTE_STATIC_EVAL(data.oYaba_augmented[joint_j].noalias()) -=
            crosscoupling_ji_Jcols_Dinv
            * crosscoupling_ji_Jcols.transpose(); // Warning: UDinv() is actually edge_ij * J, U()
                                                  // is actually edge_ij * J_cols * Dinv
                                                  //          data.of[joint_j].toVector().noalias()
                                                  //          += crosscoupling_ij * a_tmp;

          const Matrix6 crosscoupling_ji_oL = crosscoupling_ji * oL;
          if (joint_j == parent)
          {
            data.oYaba_augmented[parent].noalias() +=
              crosscoupling_ji_oL + crosscoupling_ji_oL.transpose();
          }
          else
          {
            assert(
              joint_cross_coupling.exists(JointPair(joint_j, parent))
              || joint_cross_coupling.exists(JointPair(parent, joint_j)));

            // In this particular case, the pair (joint_j,parent) might not exist, but (parent,
            // joint_j) will
            if (joint_cross_coupling.exists(JointPair(joint_j, parent)))
            {
              joint_cross_coupling.get({joint_j, parent}).noalias() += crosscoupling_ji_oL;
            }
            else
            {
              joint_cross_coupling.get({parent, joint_j}).noalias() +=
                crosscoupling_ji_oL.transpose();
            }
          }

          for (size_t k = j + 1; k < joint_neighbours.size(); ++k)
          {
            const JointIndex joint_k = joint_neighbours[k];

            assert(joint_cross_coupling.exists(JointPair(joint_k, joint_i)));
            auto & crosscoupling_ki = joint_cross_coupling.get(JointPair(joint_k, joint_i));

            assert(projected_joint_cross_coupling.exists(JointPair(joint_k, joint_i)));
            auto & crosscoupling_ki_Jcols =
              projected_joint_cross_coupling[JointPair(joint_k, joint_i)];

            PROMOTE_STATIC_EVAL(crosscoupling_ki_Jcols.noalias()) = crosscoupling_ki * Jcols;

            assert(joint_j != joint_k && "Must never happen!");

            assert(joint_cross_coupling.exists(JointPair(joint_j, joint_k)));
            auto & crosscoupling_jk = joint_cross_coupling.get(JointPair(joint_j, joint_k));
            DO_NOT_PROMOTE_STATIC_EVAL(crosscoupling_jk.noalias()) -=
              crosscoupling_ji_Jcols_Dinv * crosscoupling_ki_Jcols.transpose();
          }
        }
      }
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplApplyOnTheRightBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      const JointIndex joint_i = jmodel.id();
      const JointIndex parent = model.parents[joint_i];

      // Compare to ABA, the sign of f[joint_i] is reversed
      jmodel.jointVelocitySelector(custom_data.u) += jdata.S().transpose() * custom_data.f[joint_i];

      if (parent > 0)
      {
        auto & pa = custom_data.f[joint_i];
        // Compare to ABA, the sign of f[joint_i] is reversed
        DO_NOT_PROMOTE_STATIC_EVAL(pa.toVector().noalias()) -=
          jdata.UDinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.f[parent] += data.liMi[joint_i].act(pa);
      }
    }
  };

  template<typename DelassusOperator>
  struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass
  : public fusion::JointUnaryVisitorBase<
      DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      typedef typename Model::JointIndex JointIndex;

      const JointIndex joint_i = jmodel.id();
      const JointIndex parent = model.parents[joint_i];

      //      typename JointData::TangentVector_t ddq_joint;
      auto ddq_joint = jmodel.jointVelocitySelector(custom_data.ddq);
      if (parent > 0)
      {
        custom_data.a[joint_i] += data.liMi[joint_i].actInv(custom_data.a[parent]);
        PROMOTE_STATIC_EVAL(ddq_joint.noalias()) =
          jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u);
        PROMOTE_STATIC_EVAL(ddq_joint.noalias()) -=
          jdata.UDinv().transpose() * custom_data.a[joint_i].toVector();
        custom_data.a[joint_i] += jdata.S() * ddq_joint;
      }
      else
      {
        PROMOTE_STATIC_EVAL(ddq_joint.noalias()) =
          jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.a[joint_i] = jdata.S() * ddq_joint;
      }
    }

  }; // struct DelassusOperatorRigidBodySystemsTplApplyOnTheRightForwardPass

  template<typename DelassusOperator>
  struct AugmentedMassMatrixOperatorSolveInPlaceBackwardPass
  : public fusion::JointUnaryVisitorBase<
      AugmentedMassMatrixOperatorSolveInPlaceBackwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      typedef typename Model::Scalar Scalar;
      typedef typename Data::Force Force;
      typedef std::pair<JointIndex, JointIndex> JointPair;

      const auto & neighbours = data.joint_neighbours;
      // auto & joint_cross_coupling = data.joint_cross_coupling;
      const auto & projected_joint_cross_coupling = data.projected_joint_cross_coupling;

      const JointIndex joint_i = jmodel.id();
      const JointIndex parent = model.parents[joint_i];
      const auto & joint_neighbours = neighbours[joint_i];

      const auto Jcols = jmodel.jointCols(data.J);

      Force & ofi = custom_data.of_augmented[joint_i];

      // Compare to ABA, the sign of ofi is reversed
      PROMOTE_STATIC_EVAL(jmodel.jointVelocitySelector(custom_data.u).noalias()) +=
        Jcols.transpose() * ofi.toVector();

      if (joint_neighbours.size())
      {
        using VectorNV = typename std::remove_reference<typename JointData::TangentVector_t>::type;
        typedef Eigen::Map<VectorNV, EIGEN_DEFAULT_ALIGN_BYTES> MapVectorNV;
        MapVectorNV res = MapVectorNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, jmodel.nv(), 1));
        DO_NOT_PROMOTE_STATIC_EVAL(res.noalias()) =
          (jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u));

        // const Vector6 Ji_res = Jcols * res;

        for (JointIndex joint_j : joint_neighbours)
        {
          // const Matrix6 & crosscoupling_ji =
          //   (i > joint_j) ? joint_cross_coupling.get(JointPair(joint_j, joint_i))
          //                  : joint_cross_coupling.get(JointPair(i, joint_j)).transpose();

          assert(projected_joint_cross_coupling.exists(JointPair(joint_j, joint_i)));
          const auto & projected_crosscoupling_ji_Jcols =
            projected_joint_cross_coupling[JointPair(joint_j, joint_i)];

          Force & ofj = custom_data.of_augmented[joint_j];
          // Compare to ABA, the sign of ofj is reversed
          DO_NOT_PROMOTE_STATIC_EVAL(ofj.toVector().noalias()) -=
            projected_crosscoupling_ji_Jcols * res;
        }
      }

      if (parent > 0)
      {
        // Compare to ABA, the sign of ofi is reversed
        DO_NOT_PROMOTE_STATIC_EVAL(ofi.toVector().noalias()) -=
          jdata.UDinv() * jmodel.jointVelocitySelector(custom_data.u);
        custom_data.of_augmented[parent] += ofi;
      }
    }
  };

  template<typename DelassusOperator>
  struct AugmentedMassMatrixOperatorSolveInPlaceForwardPass
  : public fusion::JointUnaryVisitorBase<
      AugmentedMassMatrixOperatorSolveInPlaceForwardPass<DelassusOperator>>
  {
    typedef typename DelassusOperator::Model Model;
    typedef typename DelassusOperator::Data Data;
    typedef typename DelassusOperator::CustomData CustomData;
    typedef std::pair<JointIndex, JointIndex> JointPair;

    typedef boost::fusion::vector<const Model &, const Data &, CustomData &> ArgsType;

    template<typename JointModel>
    static void algo(
      const pinocchio::JointModelBase<JointModel> & jmodel,
      const pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
      const Model & model,
      const Data & data,
      CustomData & custom_data)
    {
      typedef typename Model::Scalar Scalar;
      typedef typename Model::JointIndex JointIndex;

      const auto J_cols = jmodel.jointCols(data.J);

      const JointIndex joint_i = jmodel.id();
      const JointIndex parent = model.parents[joint_i];
      const auto & joint_neighbours = data.joint_neighbours[joint_i];
      const auto & projected_joint_cross_coupling = data.projected_joint_cross_coupling;

      auto & oai = custom_data.oa_augmented[joint_i];
      oai = custom_data.oa_augmented[parent];

      if (joint_neighbours.size())
      {
        using VectorNV = typename std::remove_reference<typename JointData::TangentVector_t>::type;
        typedef Eigen::Map<VectorNV, EIGEN_DEFAULT_ALIGN_BYTES> MapVectorNV;
        MapVectorNV projected_coupling_forces =
          MapVectorNV(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, jmodel.nv(), 1));
        projected_coupling_forces.setZero();
        // Force coupling_forces = Force::Zero();

        for (const JointIndex joint_j : joint_neighbours)
        {
          // const Matrix6 & crosscoupling_ij =
          //   (i > joint_j) ? data.joint_cross_coupling.get(JointPair(joint_j,
          //   joint_i)).transpose()
          //                  : data.joint_cross_coupling.get(JointPair(i, joint_j));

          assert(projected_joint_cross_coupling.exists(JointPair(joint_j, joint_i)));
          const auto & projected_crosscoupling_ji_Jcols =
            projected_joint_cross_coupling[JointPair(joint_j, joint_i)];

          const auto & oaj = custom_data.oa_augmented[joint_j];
          // coupling_forces.toVector().noalias() += crosscoupling_ij * oaj.toVector();
          DO_NOT_PROMOTE_STATIC_EVAL(projected_coupling_forces.noalias()) +=
            projected_crosscoupling_ji_Jcols.transpose() * oaj.toVector();
        }

        jmodel.jointVelocitySelector(custom_data.u).noalias() -=
          // J_cols.transpose() * coupling_forces.toVector();
          projected_coupling_forces;
      }

      auto ddq_segment = jmodel.jointVelocitySelector(custom_data.ddq);
      PROMOTE_STATIC_EVAL(ddq_segment.noalias()) =
        jdata.Dinv() * jmodel.jointVelocitySelector(custom_data.u);
      PROMOTE_STATIC_EVAL(ddq_segment.noalias()) -= jdata.UDinv().transpose() * oai.toVector();
      DO_NOT_PROMOTE_STATIC_EVAL(oai.toVector().noalias()) += J_cols * ddq_segment;
    }
  };

} // namespace pinocchio

#undef PROMOTE_STATIC_EVAL
#undef DO_NOT_PROMOTE_STATIC_EVAL

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_visitors_hxx__
