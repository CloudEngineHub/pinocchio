//
// Copyright (c) 2016-2020 CNRS INRIA
//

#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include <Eigen/Cholesky>
  #include <Eigen/SparseCholesky>

  #include "pinocchio/context.hxx"
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

  template<typename ContactCholeskyDecomposition>
  struct DelassusCholeskyExpressionTpl;
  typedef DelassusCholeskyExpressionTpl<ContactCholeskyDecomposition> DelassusCholeskyExpression;

  template<typename Scalar, int Options>
  struct RigidConstraintModelTpl;
  template<typename Scalar, int Options>
  struct RigidConstraintDataTpl;

  typedef RigidConstraintModelTpl<context::Scalar, context::Options> RigidConstraintModel;
  typedef RigidConstraintDataTpl<context::Scalar, context::Options> RigidConstraintData;

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
} // namespace pinocchio
