//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_utils_hpp__
#define __pinocchio_algorithm_constraints_utils_hpp__

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include "pinocchio/utils/std-vector.hpp"
#include "pinocchio/utils/reference.hpp"

namespace pinocchio
{

  /**
   * @brief Create a vector of ConstraintData from a vector of ConstraintModel.
   *
   * This function iterates over a vector of constraint models and calls
   * `createData()` on each model to generate the associated constraint data.
   * The resulting data is collected into a vector that uses the same allocator
   * as the input vector.
   *
   * @tparam ConstraintModel The type of constraint model. It must define a
   *         nested type `ConstraintData` and have a method `ConstraintData createData() const`.
   * @tparam ConstraintModelAllocator The allocator used for the input vector of constraint models.
   *
   * @param constraint_models A vector of constraint models used to generate the constraint data.
   *
   * @return A vector of `ConstraintModel::ConstraintData` using the same allocator
   *         as the input vector type. Each element corresponds to the result of
   *         `createData()` called on the respective constraint model.
   *
   * @note The return type is computed using `internal::std_vector_with_same_allocator`,
   *       which ensures that the output vector has a compatible allocator with the input vector.
   */
  template<typename ConstraintModel, class ConstraintModelAllocator>
  typename internal::template std_vector_with_same_allocator<
    std::vector<ConstraintModel, ConstraintModelAllocator>>::
    template type<typename ConstraintModel::ConstraintData>
    createData(const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models)
  {
    typedef typename internal::template std_vector_with_same_allocator<
      std::vector<ConstraintModel, ConstraintModelAllocator>>::
      template type<typename ConstraintModel::ConstraintData>
        ReturnType;

    ReturnType constraint_datas;
    constraint_datas.reserve(constraint_models.size());

    for (const auto & cm : constraint_models)
      constraint_datas.push_back(cm.createData());

    return constraint_datas;
  }

  /**
   * @brief Compute the total size of a set of constraint models for the selection (Current and
   * Maximal).
   *
   * This function iterates through a list of constraint models and accumulates the number
   * of selected constraint sizes. For each constraint model it calls
   * `ConstraintModel::residualSize(sel)` to determine the number of selected constraints, then sums
   * these values over all constraints in the input vectors.
   *
   * @tparam ConstraintModel Type of each constraint model contained in the vector.
   * @tparam ConstraintModelAllocator Allocator type used for the vector of constraint models.
   *
   * @param[in] constraint_models Vector of constraint model objects.
   * @param[in] sel Either CURRENT or MAXIMAL.
   *
   * @return The total active size (dimension) obtained by summing the active sizes
   *         of all individual constraint models.
   *
   * @sa ConstraintModelTpl::residualSize
   */
  template<
    typename ConstraintModel,
    class ConstraintModelAllocator,
    ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
  Eigen::Index residualSize(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    ConstraintSelectionTag<Sel> sel = CurrentSelection())
  {
    Eigen::Index active_size = 0;
    for (std::size_t i = 0; i < constraint_models.size(); ++i)
    {
      const auto & cmodel = helper::get_ref(constraint_models[i]);
      active_size += cmodel.residualSize(sel);
    }

    return active_size;
  }

  /**
   * @brief Compute the total size of a set of constraint models for the selection (Current and
   * Maximal).
   *
   * This function iterates through a list of constraint models and accumulates the number
   * of selected constraint sizes. For each constraint model it calls
   * `ConstraintModel::symmetricConeResidualSize(sel)` to determine the number of selected
   * constraints, then sums these values over all constraints in the input vectors.
   *
   * @tparam ConstraintModel Type of each constraint model contained in the vector.
   * @tparam ConstraintModelAllocator Allocator type used for the vector of constraint models.
   *
   * @param[in] constraint_models Vector of constraint model objects.
   * @param[in] sel Either CURRENT or MAXIMAL.
   *
   * @return The total active size (dimension) obtained by summing the active sizes
   *         of all individual constraint models.
   *
   * @sa ConstraintModelTpl::symmetricConeResidualSize
   */
  template<
    typename ConstraintModel,
    class ConstraintModelAllocator,
    ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
  Eigen::Index symmetricConeResidualSize(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    ConstraintSelectionTag<Sel> sel = CurrentSelection())
  {
    Eigen::Index active_size = 0;
    for (std::size_t i = 0; i < constraint_models.size(); ++i)
    {
      const auto & cmodel = helper::get_ref(constraint_models[i]);
      active_size += cmodel.symmetricConeResidualSize(sel);
    }

    return active_size;
  }

  /**
   * @brief Compute the total size of a set of constraint models for the selection (Current and
   * Maximal).
   *
   * This function iterates through a list of constraint models and accumulates the number
   * of selected constraint sizes. For each constraint model it calls
   * `ConstraintModel::symmetricConeResidualScalingSize(sel)` to determine the number of selected
   * constraints, then sums these values over all constraints in the input vectors.
   *
   * @tparam ConstraintModel Type of each constraint model contained in the vector.
   * @tparam ConstraintModelAllocator Allocator type used for the vector of constraint models.
   *
   * @param[in] constraint_models Vector of constraint model objects.
   * @param[in] sel Either CURRENT or MAXIMAL.
   *
   * @return The total active size (dimension) obtained by summing the active sizes
   *         of all individual constraint models.
   *
   * @sa ConstraintModelTpl::symmetricConeResidualScalingSize
   */
  template<
    typename ConstraintModel,
    class ConstraintModelAllocator,
    ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
  Eigen::Index symmetricConeResidualScalingSize(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    ConstraintSelectionTag<Sel> sel = CurrentSelection())
  {
    Eigen::Index active_size = 0;
    for (std::size_t i = 0; i < constraint_models.size(); ++i)
    {
      const auto & cmodel = helper::get_ref(constraint_models[i]);
      active_size += cmodel.symmetricConeResidualScalingSize(sel);
    }

    return active_size;
  }

  template<
    typename ConstraintModel,
    class ConstraintModelAllocator,
    typename ComplianceVector,
    ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
  void setConstraintCompliance(
    std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const Eigen::MatrixBase<ComplianceVector> & compliance,
    ConstraintSelectionTag<Sel> sel = CurrentSelection())
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(ComplianceVector);

    Eigen::Index constraint_index = 0;

    assert(compliance.size() == residualSize(constraint_models, sel));

    for (std::size_t i = 0; i < constraint_models.size(); i++)
    {
      auto & cmodel = helper::get_ref(constraint_models[i]);
      const auto csize = cmodel.residualSize(sel);
      cmodel.setCompliance(compliance.segment(constraint_index, csize), sel);
      constraint_index += csize;
    }
  }

  template<
    typename ConstraintModel,
    class ConstraintModelAllocator,
    typename ComplianceVector,
    ConstraintSelectionType Sel = ConstraintSelectionType::CURRENT>
  void retrieveConstraintCompliance(
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const Eigen::MatrixBase<ComplianceVector> & compliance_,
    ConstraintSelectionTag<Sel> sel = CurrentSelection())
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(ComplianceVector);

    Eigen::Index constraint_index = 0;
    auto & compliance = compliance_.const_cast_derived();

    assert(compliance.size() == residualSize(constraint_models, sel));

    for (std::size_t i = 0; i < constraint_models.size(); i++)
    {
      const auto & cmodel = helper::get_ref(constraint_models[i]);
      const auto csize = cmodel.residualSize(sel);
      cmodel.retrieveCompliance(compliance.segment(constraint_index, csize), sel);
      constraint_index += csize;
    }
  }

  /**
   * @brief Evaluate all the constraint models given a Pinocchio model and data.
   *
   * This function iterates through the provided list of constraint models
   * and calls their respective `calc` methods using the given model and data.
   * The computation results for each constraint model are stored in the
   * corresponding constraint data objects.
   *
   * @tparam Scalar Type of the scalar values (typically double or float).
   * @tparam Options Template options passed to Pinocchio's ModelTpl and DataTpl.
   * @tparam JointCollectionTpl The joint collection type defining the set of joints supported by
   * the model.
   * @tparam ConstraintModel Type of each constraint model contained in the @p constraint_models
   * vector.
   * @tparam ConstraintModelAllocator Allocator type for the @p constraint_models vector.
   * @tparam ConstraintData Type of each constraint model contained in the @p constraint_models
   * vector.
   * @tparam ConstraintDataAllocator Allocator type for the @p constraint_datas vector.
   *
   * @param[in] model The Pinocchio model structure.
   * @param[in] data The Pinocchio data structure associated with the model.
   * @param[in] constraint_models A vector of constraint model objects describing the constraints to
   * be evaluated.
   * @param[out] constraint_datas A vector of constraint data objects where each corresponding
   * constraint model’s results are stored.
   *
   * @note The size of @p constraint_models and @p constraint_datas must match.
   * @warning This function assumes that each constraint model and its corresponding data object
   * refer to the same type of constraint.
   *
   * @sa ConstraintModelTpl::calc, ConstraintDataTpl
   */
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModel,
    class ConstraintModelAllocator,
    typename ConstraintData,
    typename ConstraintDataAllocator>
  void calc(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    for (size_t k = 0; k < constraint_models.size(); ++k)
    {
      const auto & cmodel = helper::get_ref(constraint_models[k]);
      auto & cdata = helper::get_ref(constraint_datas[k]);

      cmodel.calc(model, data, cdata);
    }
  }

  ///
  /// \brief Computes the kinematic Jacobian associatied to a given constraint model.
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  /// This function also assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_model Constraint model.
  /// \param[in] constraint_data Constraint data.
  /// \param[out] J A reference on the Jacobian matrix where the results will be stored in (dim 6 x
  /// model.nv). You must fill J with zero elements, e.g. J.fill(0.).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModelDerived,
    typename ConstraintDataDerived,
    typename Matrix6Like>
  void getConstraintJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const ConstraintModelBase<ConstraintModelDerived> & constraint_model,
    const ConstraintDataBase<ConstraintDataDerived> & constraint_data,
    const Eigen::MatrixBase<Matrix6Like> & J);

  ///
  /// \brief Computes the kinematic Jacobian associatied to a given set of constraint models.
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  /// This function also assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \param[out] J A reference on the Jacobian matrix where the results will be stored in (dim nc x
  /// model.nv). You must fill J with zero elements, e.g. J.fill(0.).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename DynamicMatrixLike>
  void getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<DynamicMatrixLike> & J);

  ///
  /// \brief Computes the kinematic Jacobian associatied to a given set of constraint models.
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  /// This function also assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \return A Jacobian matrix where the results will be stored in (dim nc x model.nv).
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  typename DataTpl<Scalar, Options, JointCollectionTpl>::MatrixXs getConstraintsJacobian(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

  ///
  /// \brief Evaluate the operation res = J * rhs
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  /// This function also assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \param[in] rhs Right-hand side term.
  /// \param[out] res Results.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType,
    AssignmentOperatorType op = SETTO>
  void evalConstraintJacobianMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res,
    AssignmentOperatorTag<op> aot = SetTo());

  ///
  /// \brief Evaluate the operation res = J.T * rhs
  ///
  /// \remarks This function assumes that the a computeJointJacobians has been called first or any
  /// algorithms that computes data.J and data.oMi.
  /// This function also assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint data.
  /// \param[in] rhs Right-hand side term.
  /// \param[out] res Results.
  ///
  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename RhsMatrixType,
    typename ResultMatrixType,
    AssignmentOperatorType op = SETTO>
  void evalConstraintJacobianTransposeMatrixProduct(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<RhsMatrixType> & rhs,
    const Eigen::MatrixBase<ResultMatrixType> & res,
    AssignmentOperatorTag<op> aot = SetTo());

  ///
  /// \brief Maps the constraint forces expressed in the constraint space to joint forces expressed
  /// in the local frame.
  ///
  /// \remarks This function assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint datas.
  /// \param[in] constraint_forces Matrix or vector containing the constraint forces.
  /// \param[out] joint_forces Vector of  joint forces (dimension model.njoints).
  ///
  template<
    typename Scalar,
    int Options,
    int ForceOptions,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename ForceMatrix,
    class ForceAllocator,
    ReferenceFrame rf>
  void mapConstraintForcesToJointForces(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<ForceMatrix> & constraint_forces,
    std::vector<ForceTpl<Scalar, ForceOptions>, ForceAllocator> & joint_forces,
    ReferenceFrameTag<rf> reference_frame);

  ///
  /// \brief Maps the joint motions expressed in the joint space local frame to the constraint
  /// motions.
  ///
  /// \remarks This function assumes that the constrained datas are up-to-date.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] constraint_models Vector of constraint models.
  /// \param[in] constraint_datas Vector of constraint datas.
  /// \param[in] joint_motions Vector of  joint motions (dimension model.njoints).
  /// \param[out] constraint_motions Resulting matrix or vector containing the constraint motions.
  ///
  template<
    typename Scalar,
    int Options,
    int MotionOptions,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    class MotionAllocator,
    typename MotionMatrix,
    ReferenceFrame rf>
  void mapJointMotionsToConstraintMotions(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const std::vector<MotionTpl<Scalar, MotionOptions>, MotionAllocator> & joint_motions,
    const Eigen::MatrixBase<MotionMatrix> & constraint_motions,
    ReferenceFrameTag<rf> reference_frame);
} // namespace pinocchio

#include "pinocchio/algorithm/constraints/utils.hxx"

#if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #include "pinocchio/algorithm/constraints/utils.txx"
#endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#endif // __pinocchio_algorithm_constraints_utils_hpp__
