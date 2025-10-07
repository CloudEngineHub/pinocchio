//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_constraints_utils_hpp__
#define __pinocchio_algorithm_constraints_utils_hpp__

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

  template<
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    typename ConstraintModel,
    class ConstraintModelAllocator,
    typename ConstraintData,
    class ConstraintDataAllocator>
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

} // namespace pinocchio

#endif // __pinocchio_algorithm_constraints_utils_hpp__
