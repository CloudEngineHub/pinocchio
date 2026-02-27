//
// Copyright (c) 2016-2020 CNRS INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include <vector>

  #include <Eigen/Cholesky>
  #include <Eigen/SparseCholesky>

  #include "pinocchio/context.hpp"
  #include "pinocchio/multibody/fwd.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<typename Scalar>
  struct ProximalSettingsTpl;
  typedef ProximalSettingsTpl<context::Scalar> ProximalSettings;

  template<typename Scalar, int Options>
  struct ContactCholeskyDecompositionTpl;
  typedef ContactCholeskyDecompositionTpl<context::Scalar, context::Options>
    ContactCholeskyDecomposition;

  template<typename Scalar, int Options>
  struct RigidConstraintModelTpl;
  template<typename Scalar, int Options>
  struct RigidConstraintDataTpl;

  typedef RigidConstraintModelTpl<context::Scalar, context::Options> RigidConstraintModel;
  typedef RigidConstraintDataTpl<context::Scalar, context::Options> RigidConstraintData;
  typedef std::vector<RigidConstraintModel> RigidConstraintModelVector;
  typedef std::vector<RigidConstraintData> RigidConstraintDataVector;

  template<typename DelassusOperatorDerived>
  struct DelassusOperatorBase;

  template<typename DelassusOperatorDerived, typename MatrixDerived>
  struct DelassusOperatorApplyOnTheRightReturnType;

  template<
    typename Scalar,
    int Options = 0,
    template<typename, auto...> class CholeskyDecompositionTpl = Eigen::LLT>
  struct DelassusOperatorDenseTpl;
  typedef DelassusOperatorDenseTpl<context::Scalar, context::Options> DelassusOperatorDense;

  template<
    typename Scalar,
    int Options = 0,
    class SparseCholeskyDecomposition = Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar, Options>>>
  struct DelassusOperatorSparseTpl;
  typedef DelassusOperatorSparseTpl<context::Scalar, context::Options> DelassusOperatorSparse;

  template<
    typename _Scalar,
    int _Options,
    template<typename, int> class _JointCollectionTpl,
    class _ConstraintModel,
    template<typename T> class Holder>
  struct DelassusOperatorRigidBodySystemsTpl;

  template<typename DelassusOperator, typename PreconditionerType>
  struct DelassusOperatorPreconditionedTpl;

  template<typename ContactCholeskyDecomposition>
  struct DelassusCholeskyExpressionTpl;
  typedef DelassusCholeskyExpressionTpl<ContactCholeskyDecomposition> DelassusCholeskyExpression;

  template<class... D>
  struct AlgorithmCheckerList;

  struct ParentChecker;
  struct CRBAChecker;
  struct ABAChecker;

  AlgorithmCheckerList<ParentChecker, CRBAChecker, ABAChecker> makeDefaultCheckerList();

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  bool checkData(
    const ModelTpl<Scalar, Options, JointCollectionTpl> & model,
    const DataTpl<Scalar, Options, JointCollectionTpl> & data);

  template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options>
  neutral(const ModelTpl<Scalar, Options, JointCollectionTpl> & model);
} // namespace pinocchio
