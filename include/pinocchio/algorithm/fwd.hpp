//
// Copyright (c) 2020-2025 INRIA
//

#ifndef __pinocchio_algorithm_fwd_hpp__
#define __pinocchio_algorithm_fwd_hpp__

#include "pinocchio/fwd.hpp"

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

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModel,
    template<typename T> class Holder = std::reference_wrapper>
  struct DelassusOperatorRigidBodySystemsTpl;
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_fwd_hpp__
