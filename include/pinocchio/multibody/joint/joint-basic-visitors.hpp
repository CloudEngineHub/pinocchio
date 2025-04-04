//
// Copyright (c) 2016-2020 CNRS INRIA
//

#ifndef __pinocchio_multibody_joint_basic_visitors_hpp__
#define __pinocchio_multibody_joint_basic_visitors_hpp__

#include "pinocchio/multibody/joint/fwd.hpp"

namespace pinocchio
{

  /**
   * @brief      Visit a JointModelTpl through CreateData visitor to create a JointDataTpl.
   *
   * @tparam JointCollection    Collection of Joint types.
   *
   * @param[in]  jmodel  The JointModelTpl we want to create a data for.
   *
   * @return     The created JointDataTpl
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline JointDataTpl<Scalar, Options, JointCollectionTpl>
  createData(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl and the corresponding JointDataTpl through
   * JointCalcZeroOrderVisitor to compute the joint data kinematics at order zero.
   *
   * @tparam JointCollection    Collection of Joint types.
   * @tparam ConfigVectorType   Type of the joint configuration vector.
   *
   * @param[in]  jmodel  The corresponding JointModelVariant to the JointDataVariant we want to
   * update
   * @param      jdata   The JointDataVariant we want to update
   * @param[in]  q       The full model's (in which the joint belongs to) configuration vector
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename ConfigVectorType>
  inline void calc_zero_order(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
    JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata,
    const Eigen::MatrixBase<ConfigVectorType> & q);

  /**
   * @brief      Visit a JointModelTpl and the corresponding JointDataTpl through
   * JointCalcFirstOrderVisitor to compute the joint data kinematics at order one.
   *
   * @tparam JointCollection    Collection of Joint types.
   * @tparam ConfigVectorType   Type of the joint configuration vector.
   * @tparam TangentVectorType  Type of the joint velocity vector.
   *
   * @param[in]  jmodel  The corresponding JointModelVariant to the JointDataVariant we want to
   * update
   * @param      jdata   The JointDataVariant we want to update
   * @param[in]  q       The full model's (in which the joint belongs to) configuration vector
   * @param[in]  v       The full model's (in which the joint belongs to) velocity vector
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename ConfigVectorType,
    typename TangentVectorType>
  inline void calc_first_order(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
    JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata,
    const Eigen::MatrixBase<ConfigVectorType> & q,
    const Eigen::MatrixBase<TangentVectorType> & v);

  /**
   * @brief      Visit a JointModelTpl and the corresponding JointDataTpl through
   * JointCalcFirstOrderVisitor to compute the joint data kinematics at order one.
   *
   * @tparam JointCollection    Collection of Joint types.
   * @tparam TangentVectorType  Type of the joint velocity vector.
   *
   * @param[in]  jmodel  The corresponding JointModelVariant to the JointDataVariant we want to
   * update
   * @param      jdata   The JointDataVariant we want to update
   * @param[in]  v       The full model's (in which the joint belongs to) velocity vector
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename TangentVectorType>
  inline void calc_first_order(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
    JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata,
    const Blank blank,
    const Eigen::MatrixBase<TangentVectorType> & v);

  /**
   * @brief      Visit a JointModelTpl and the corresponding JointDataTpl through
   * JointCalcAbaVisitor to.
   *
   * @tparam JointCollection    Collection of Joint types.
   * @tparam Matrix6Type        A matrix 6x6 like Eigen container.
   *
   * @param[in]         jmodel      The corresponding JointModelVariant to the JointDataVariant we
   * want to update
   * @param[in,out]     jdata        The JointDataVariant we want to update
   * @param[in]         armature  Armature related to the current joint.
   * @param[in,out]     I                 Inertia matrix of the subtree following the jmodel in the
   * kinematic chain as dense matrix   *
   * @param[in]         update_I  If I should be updated or not
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename VectorLike,
    typename Matrix6Type>
  inline void calc_aba(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
    JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata,
    const Eigen::MatrixBase<VectorLike> & armature,
    const Eigen::MatrixBase<Matrix6Type> & I,
    const bool update_I);

  /**
   * @brief      Visit a JointModelTpl through JointNvVisitor to get the dimension of
   *             the joint tangent space
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The dimension of joint tangent space
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline int nv(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointNqVisitor to get the dimension of
   *             the joint configuration space
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The dimension of joint configuration space
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline int nq(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointNvExtendVisitor to get the dimension of
   *             the joint extended tangent space
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The dimension of joint extended tangent space
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline int nvExtended(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointConfigurationLimitVisitor
   *             to get the configurations limits
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The bool with configurations limits of the joint
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline const std::vector<bool>
  hasConfigurationLimit(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointConfigurationLimitInTangentVisitor
   *             to get the configurations limits in tangent space
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The bool with configurations limits in tangent space of the joint
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline const std::vector<bool>
  hasConfigurationLimitInTangent(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointIdxQVisitor to get the index in the full model
   * configuration space corresponding to the first degree of freedom of the Joint
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The index in the full model configuration space corresponding to the first
   *             degree of freedom of jmodel
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline int idx_q(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointIdxVVisitor to get the index in the model
   * tangent space corresponding to the first joint tangent space degree
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The index in the full model tangent space corresponding to the first
   *             joint tangent space degree
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline int idx_v(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointIdvExtendedVisitor to get the index in the model
   * extended tangent space corresponding to the joint first joint extended tangent space degree
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The index in the model extended tangent space corresponding to the first
   *             joint extended tangent space degree
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline int idx_vExtended(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointIdVisitor to get the index of the joint in the
   * kinematic chain
   *
   * @param[in]  jmodel  The JointModelVariant
   *
   * @return     The index of the joint in the kinematic chain
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline JointIndex id(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl through JointSetIndexesVisitor to set
   *             the indexes of the joint in the kinematic chain
   *
   * @param[in]  jmodel  The JointModelVariant
   * @param[in]  id      The index of joint in the kinematic chain
   * @param[in]  q       The index in the full model configuration space corresponding to the first
   * degree of freedom
   * @param[in]  v       The index in the full model tangent space corresponding to the first joint
   * tangent space degree
   * @param[in]  vExtended       The index in the model extended tangent space corresponding to the
   * joint first extended tangent space degree
   *
   * @return     The index of the joint in the kinematic chain
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline void setIndexes(
    JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
    JointIndex id,
    int q,
    int v,
    int vExtended);

  /**
   * @brief      Visit a JointModelTpl through JointSetIndexesVisitor to set
   *             the indexes of the joint in the kinematic chain
   *
   * @param[in]  jmodel  The JointModelVariant
   * @param[in]  id      The index of joint in the kinematic chain
   * @param[in]  q       The index in the full model configuration space corresponding to the first
   * degree of freedom
   * @param[in]  v       The index in the full model tangent space corresponding to the first joint
   * tangent space degree
   *
   * @return     The index of the joint in the kinematic chain
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline void setIndexes(
    JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel, JointIndex id, int q, int v);

  /**
   * @brief      Visit a JointModelTpl through JointShortnameVisitor to get the shortname of the
   * derived joint model
   *
   * @param      jmodel  The JointModelVariant we want the shortname of the type held in
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline std::string shortname(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl<Scalar,...> to cast it into JointModelTpl<NewScalar,...>
   *
   * @tparam     NewScalar new scalar type of of the JointModelTpl
   *
   * @param[in]  jmodel  The joint model to cast.
   *
   * @return     A new JointModelTpl<NewScalar,...> casted from JointModelTpl<Scalar,...>.
   */
  template<
    typename NewScalar,
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl>
  typename CastType<NewScalar, JointModelTpl<Scalar, Options, JointCollectionTpl>>::type
  cast_joint(const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel);

  /**
   * @brief      Visit a JointModelTpl<Scalar,...> to compare it to JointModelDerived
   *
   * @param[in]  jmodel_generic  The generic joint model containing a variant.
   * @param[in]  jmodel  The other joint model for the comparison.
   *
   * @return     True if the two joint models are equal.
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename JointModelDerived>
  bool isEqual(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel_generic,
    const JointModelBase<JointModelDerived> & jmodel);

  /**
   * @brief      Check whether JointModelTpl<Scalar,...> has the indexes than another
   * JointModelDerived
   *
   * @param[in]  jmodel_generic  The generic joint model containing a variant.
   * @param[in]  jmodel  The other joint modelto compare with
   *
   * @return     True if the two joints have the same indexes.
   */
  template<
    typename NewScalar,
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename JointModelDerived>
  bool hasSameIndexes(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel_generic,
    const JointModelBase<JointModelDerived> & jmodel);

  //
  // Visitors on JointDatas
  //

  /**
   * @brief      Visit a JointDataVariant through JointConfigVisitor to get the joint configuration
   * vector
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The current value of the joint configuration vector
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline typename JointDataTpl<Scalar, Options, JointCollectionTpl>::ConfigVector_t
  joint_q(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataVariant through JointConfigVisitor to get the joint velocity
   * vector
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The current value of the joint velocity vector
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline typename JointDataTpl<Scalar, Options, JointCollectionTpl>::TangentVector_t
  joint_v(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataVariant through JointConstraintVisitor to get the joint constraint
   *             as a dense constraint
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The constraint dense corresponding to the joint derived constraint
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline JointMotionSubspaceTpl<Eigen::Dynamic, Scalar, Options>
  joint_motin_subspace_xd(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointTransformVisitor to get the joint internal
   * transform  (transform between the entry frame and the exit frame of the joint).
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The joint transform corresponding to the joint derived transform (sXp)
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline SE3Tpl<Scalar, Options>
  joint_transform(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointMotionVisitor to get the joint internal motion
   *             as a dense motion.
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The motion dense corresponding to the joint derived motion
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline MotionTpl<Scalar, Options>
  motion(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointBiasVisitor to get the joint bias
   *             as a dense motion.
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The motion dense corresponding to the joint derived bias
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline MotionTpl<Scalar, Options>
  bias(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointUInertiaVisitor to get the U matrix of the
   * inertia matrix decomposition.
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The U matrix of the inertia matrix decomposition
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options>
  u_inertia(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointDInvInertiaVisitor to get the D^{-1} matrix of
   * the inertia matrix decomposition.
   *
   * @param[in]  jdata  The jdata
   *
   * @return     The D^{-1} matrix of the inertia matrix decomposition
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options>
  dinv_inertia(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointUDInvInertiaVisitor to get U*D^{-1} matrix of the
   * inertia matrix decomposition.
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The U*D^{-1} matrix of the inertia matrix decomposition
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline Eigen::Matrix<Scalar, 6, Eigen::Dynamic, Options>
  udinv_inertia(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl through JointStUInertiaVisitor to get St*I*S matrix of the
   * inertia matrix decomposition.
   *
   * @param[in]  jdata  The joint data to visit.
   *
   * @return     The St*I*S matrix
   */
  template<typename Scalar, int Options, template<typename S, int O> class JointCollectionTpl>
  inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options>
  stu_inertia(const JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata);

  /**
   * @brief      Visit a JointDataTpl<Scalar,...> to compare it to another JointData
   *
   * @param[in]  jdata_generic  The generic joint data containing a variant.
   * @param[in]  jdata  The other joint data for the comparison.
   *
   * @return     True if the two joints data are equal.
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename JointDataDerived>
  bool isEqual(
    const JointDataTpl<Scalar, Options, JointCollectionTpl> & jmodel_generic,
    const JointDataBase<JointDataDerived> & jmodel);

  /**
   * @brief Apply the correct affine transform, on a joint configuration, depending on the joint
   * type.
   *
   * @tparam Scalar Type of scaling and offset scalars.
   * @tparam Options
   * @tparam JointCollectionTpl Collection of Joint types
   * @tparam ConfigVectorIn Type of the input joint configuration vector.
   * @tparam ConfigVectorOut Type of the ouptut joint configuration vector.
   * @param jmodel Joint variant to determine the correct affine transform to use.
   * @param qIn Input configuration vector
   * @param scaling scaling factor
   * @param offset Offset value
   * @param qOut Ouptut joint configuration vector
   */
  template<
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename ConfigVectorIn,
    typename ConfigVectorOut>
  void configVectorAffineTransform(
    const JointModelTpl<Scalar, Options, JointCollectionTpl> & jmodel,
    const Eigen::MatrixBase<ConfigVectorIn> & qIn,
    const Scalar & scaling,
    const Scalar & offset,
    const Eigen::MatrixBase<ConfigVectorOut> & qOut);

  /**
   * @brief Apply joint constraint (motion subspace) on a force.
   * @tparam Op AssignmentOperatorType.
   * @tparam Scalar Type of scaling and offset scalars.
   * @tparam Options.
   * @tparam JointCollectionTpl Collection of Joint types.
   * @tparam ForceType Type of the force to transform.
   * @tparam ExpressionType Eigen expression used to store the result.
   * @param jdata Joint data to extract the motion subspace (S).
   * @param F Force applied on the motion subspace.
   * @param R Eigen expression to store the result.
   */
  template<
    int Op,
    typename Scalar,
    int Options,
    template<typename S, int O> class JointCollectionTpl,
    typename ForceType,
    typename ExpressionType>
  void applyConstraintOnForceVisitor(
    JointDataTpl<Scalar, Options, JointCollectionTpl> & jdata, ForceType F, ExpressionType R);

} // namespace pinocchio

/* --- Details -------------------------------------------------------------------- */
// Included later
// #include "pinocchio/multibody/joint/joint-basic-visitors.hxx"

#endif // ifndef __pinocchio_multibody_joint_basic_visitors_hpp__
