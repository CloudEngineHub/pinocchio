//
// Copyright (c) 2016-2020 CNRS INRIA
//
 
#pragma once

#include "pinocchio/multibody/fwd.hxx"

namespace pinocchio
{

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