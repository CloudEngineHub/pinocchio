//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_algorithm_constraint_cholesky_hxx__
#define __pinocchio_algorithm_constraint_cholesky_hxx__

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/utils/reference.hpp"
#include "pinocchio/utils/size-in-bytes.hpp"

#include "pinocchio/algorithm/constraints/constraints.hpp"

#include <algorithm>
#include <cstddef>

namespace pinocchio
{

  // TODO Remove when API is stabilized
  PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
  PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS

  template<typename Scalar, int Options>
  ContactCholeskyDecompositionTpl<Scalar, Options>::ContactCholeskyDecompositionTpl()
  : D(D_storage.map())
  , Dinv(Dinv_storage.map())
  , U(U_storage.map())
  , compliance(compliance_storage.map())
  , damping(damping_storage.map())
  , delassus_block(delassus_block_storage.map())
  {
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
  ContactCholeskyDecompositionTpl<Scalar, Options>::ContactCholeskyDecompositionTpl(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    const DataTpl<S1, O1, JointCollectionTpl> & data)
  : D(D_storage.map())
  , Dinv(Dinv_storage.map())
  , U(U_storage.map())
  , compliance(compliance_storage.map())
  , damping(damping_storage.map())
  , delassus_block(delassus_block_storage.map())
  {
    typedef ConstraintModelTpl<Scalar, Options> ConstraintModel;
    typedef ConstraintDataTpl<Scalar, Options> ConstraintData;

    std::vector<ConstraintModel> empty_constraint_models;
    std::vector<ConstraintData> empty_constraint_datas;
    resize(model, data, empty_constraint_models, empty_constraint_datas);
  }

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  ContactCholeskyDecompositionTpl<Scalar, Options>::ContactCholeskyDecompositionTpl(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    const DataTpl<S1, O1, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  : D(D_storage.map())
  , Dinv(Dinv_storage.map())
  , U(U_storage.map())
  , compliance(compliance_storage.map())
  , damping(damping_storage.map())
  , delassus_block(delassus_block_storage.map())
  {
    PINOCCHIO_UNUSED_VARIABLE(data);
    resize(model, data, constraint_models, constraint_datas);
  }

  template<typename Scalar, int Options>
  ContactCholeskyDecompositionTpl<Scalar, Options>::ContactCholeskyDecompositionTpl(
    const ContactCholeskyDecompositionTpl & other)
  : D(D_storage.map())
  , Dinv(Dinv_storage.map())
  , U(U_storage.map())
  , compliance(compliance_storage.map())
  , damping(damping_storage.map())
  , delassus_block(delassus_block_storage.map())
  {
    *this = other;
  }

  template<typename Scalar, int Options>
  ContactCholeskyDecompositionTpl<Scalar, Options> &
  ContactCholeskyDecompositionTpl<Scalar, Options>::operator=(
    const ContactCholeskyDecompositionTpl & other)
  {
    parents_fromRow = other.parents_fromRow;
    nv_subtree_fromRow = other.nv_subtree_fromRow;
    nv = other.nv;

    rowise_sparsity_pattern = other.rowise_sparsity_pattern;

    D_storage = other.D_storage;
    Dinv_storage = other.Dinv_storage;
    U_storage = other.U_storage;
    compliance_storage = other.compliance_storage;
    damping_storage = other.damping_storage;

    return *this;
  }

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  PINOCCHIO_DEPRECATED void ContactCholeskyDecompositionTpl<Scalar, Options>::allocate(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    const DataTpl<S1, O1, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    resize(model, data, constraint_models, constraint_datas);
  }

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::resize(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    const DataTpl<S1, O1, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas)
  {
    assert(
      constraint_models.size() == constraint_datas.size()
      && "Both std::vector should be of equal size.");
    assert(model.check(MimicChecker()) && "Function does not support mimic joints");

    nv = model.nv;
    const auto total_constraint_size = residualSize(constraint_models, constraint_datas);

    const Eigen::Index total_size = nv + total_constraint_size;

    // Compute first parents_fromRow for all the joints.
    // This code is very similar to the code of Data::computeParents_fromRow,
    // but shifted with a value corresponding to the number of constraints.
    parents_fromRow.resize(total_size);
    parents_fromRow.fill(-1);

    nv_subtree_fromRow.resize(total_size);
    //      nv_subtree_fromRow.fill(0);

    // Fill nv_subtree_fromRow for model
    for (JointIndex joint_id = 1; joint_id < (JointIndex)(model.njoints); joint_id++)
    {
      const JointIndex parent_id = model.parents[joint_id];

      const auto & joint = model.joints[joint_id];
      const auto & parent_joint = model.joints[parent_id];
      const int nvj = joint.nv();
      const int idx_vj = joint.idx_v();

      if (parent_id > 0)
        parents_fromRow[idx_vj + total_constraint_size] =
          parent_joint.idx_v() + parent_joint.nv() - 1 + total_constraint_size;
      else
        parents_fromRow[idx_vj + total_constraint_size] = -1;

      const JointIndex last_child =
        model.subtrees[joint_id].size() > 0 ? model.subtrees[joint_id].back() : JointIndex(0);
      nv_subtree_fromRow[idx_vj + total_constraint_size] =
        model.joints[last_child].idx_v() + model.joints[last_child].nv() - idx_vj;

      for (int row = 1; row < nvj; ++row)
      {
        parents_fromRow[idx_vj + total_constraint_size + row] =
          idx_vj + row - 1 + total_constraint_size;
        nv_subtree_fromRow[idx_vj + total_constraint_size + row] =
          nv_subtree_fromRow[idx_vj + total_constraint_size] - row;
      }
    }

    Eigen::Index row_id = 0;
    for (std::size_t i = 0; i < constraint_models.size(); i++)
    {
      const auto & cmodel = helper::get_ref(constraint_models[i]);
      const auto & cdata = helper::get_ref(constraint_datas[i]);
      for (Eigen::Index k = 0; k < cmodel.residualSize(cdata); ++k, row_id++)
      {
        const auto & row_active_indexes = cmodel.getRowIndexes(model, data, cdata, k);
        nv_subtree_fromRow[row_id] =
          total_constraint_size - row_id + 1
          + (row_active_indexes.size() > 0 ? row_active_indexes.back() : 0);
      }
    }
    assert(row_id == total_constraint_size);

    // Fill the sparsity pattern for each Row of the Cholesky decomposition (matrix U)
    /*
          static const Slice default_slice_value(1,1);
          static const SliceVector default_slice_vector(1,default_slice_value);

          rowise_sparsity_pattern.clear();
          rowise_sparsity_pattern.resize((size_t)total_constraint_size,default_slice_vector);
          row_id = 0; size_t constraint_id = 0;
          for(typename RigidConstraintModelVector::const_iterator it = constraint_models.begin();
              it != constraint_models.end();
              ++it, ++constraint_id)
          {
            const RigidConstraintModel & cmodel = *it;
            const BooleanVector & joint1_indexes_ee = cmodel.colwise_joint1_sparsity;
            const Eigen::Index contact_dim = cmodel.size();

            for(Eigen::Index k = 0; k < contact_dim; ++k)
            {
              SliceVector & slice_vector = rowise_sparsity_pattern[(size_t)row_id];
              slice_vector.clear();
              slice_vector.push_back(Slice(row_id,total_constraint_size-row_id));

              bool previous_index_was_true = true;
              for(Eigen::Index joint1_indexes_constraint_id = total_constraint_size;
                  joint1_indexes_constraint_id < total_size;
                  ++joint1_indexes_constraint_id)
              {
                if(joint1_indexes_ee[joint1_indexes_constraint_id])
                {
                  if(previous_index_was_true) // no discontinuity
                    slice_vector.back().size++;
                  else // discontinuity; need to create a new slice
                  {
                    const Slice new_slice(joint1_indexes_constraint_id,1);
                    slice_vector.push_back(new_slice);
                  }
                }

                previous_index_was_true = joint1_indexes_ee[joint1_indexes_constraint_id];
              }

              row_id++;
            }
          }
     */

    // Allocate Eigen memory if needed
    compliance_storage.resize(total_constraint_size);
    compliance.setZero();
    damping_storage.resize(total_constraint_size);
    damping.setZero();

    D_storage.resize(total_size);
    Dinv_storage.resize(total_size);
    U_storage.resize(total_size, total_size);
    delassus_block_storage.resize(total_constraint_size, total_constraint_size);
    U.setIdentity();
  }

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator,
    typename VectorLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::compute(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    DataTpl<S1, O1, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const Eigen::MatrixBase<VectorLike> & mus)
  {
    assert(model.check(data) && "data is not consistent with model.");
    assert(model.check(MimicChecker()) && "Function does not support mimic joints");

    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      constraint_models.size() == constraint_datas.size(),
      "The number of constraints between constraint_models and constraint_datas vectors is "
      "different.");
    PINOCCHIO_ONLY_USED_FOR_DEBUG(model);

    const Eigen::Index total_constraint_size = constraintDim();

    const auto & M = data.M;

    const size_t num_constraints = constraint_models.size();

    // Fill the mass matrix part
    D.tail(model.nv) = M.diagonal();
    U.bottomRightCorner(model.nv, model.nv).template triangularView<Eigen::StrictlyUpper>() =
      M.template triangularView<Eigen::StrictlyUpper>();

    // Constraint filling
    Eigen::Index current_row = 0;
    U.topRightCorner(total_constraint_size, model.nv).setZero();
    for (size_t constraint_id = 0; constraint_id < num_constraints; ++constraint_id)
    {
      const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
      const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);

      const Eigen::Index constraint_size = cmodel.residualSize(cdata);
      auto U_block = U.block(current_row, total_constraint_size, constraint_size, model.nv);
      cmodel.jacobian(model, data, cdata, U_block);
      current_row += constraint_size;
    }

    // Cholesky decomposition
    for (Eigen::Index j = nv - 1; j >= 0; --j)
    {
      // Classic Cholesky decomposition related to the mass matrix
      const Eigen::Index jj = total_constraint_size + j; // shifted index
      const Eigen::Index NVT = nv_subtree_fromRow[jj] - 1;

      typedef Eigen::Map<Vector, EIGEN_DEFAULT_ALIGN_BYTES> MapVector;
      MapVector DUt_partial = MapVector(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, NVT, 1));

      if (NVT)
        DUt_partial.noalias() =
          U.row(jj).segment(jj + 1, NVT).transpose().cwiseProduct(D.segment(jj + 1, NVT));

      D[jj] -= U.row(jj).segment(jj + 1, NVT).dot(DUt_partial);
      assert(
        check_expression_if_real<Scalar>(D[jj] != Scalar(0))
        && "The diagonal element is equal to zero.");
      Dinv[jj] = Scalar(1) / D[jj];

      for (Eigen::Index _ii = parents_fromRow[jj]; _ii >= total_constraint_size;
           _ii = parents_fromRow[_ii])
      {
        U(_ii, jj) -= U.row(_ii).segment(jj + 1, NVT).dot(DUt_partial);
        U(_ii, jj) *= Dinv[jj];
      }

      // Constraint part
      Eigen::Index current_row = total_constraint_size - 1;
      for (size_t index = 0; index < num_constraints; ++index)
      {
        const size_t constraint_id = num_constraints - 1 - index;
        const auto & cmodel = helper::get_ref(constraint_models[constraint_id]);
        const auto & cdata = helper::get_ref(constraint_datas[constraint_id]);
        const Eigen::Index constraint_size = cmodel.residualSize(cdata);

        for (Eigen::Index constraint_row_id = constraint_size - 1; constraint_row_id >= 0;
             --constraint_row_id, --current_row)
        {
          const auto & colwise_sparsity =
            cmodel.getRowSparsityPattern(model, data, cdata, constraint_row_id);
          if (colwise_sparsity[j])
          {
            U(current_row, jj) -= U.row(current_row).segment(jj + 1, NVT).dot(DUt_partial);
            U(current_row, jj) *= Dinv[jj];
          }
        }
      }
    }

    // Setting physical compliance
    retrieveCompliance(constraint_models, constraint_datas, compliance);

    // Setting numerical damping
    {
      computeDelassusBlock();
      updateDamping(mus);
    }
  }

  template<typename Scalar, int Options>
  template<typename VectorLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateCompliance(
    const Eigen::MatrixBase<VectorLike> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
    compliance = vec;

    // The diagonal term of the KKT should be updated with the new compliance
    updateDamping(getDamping());
  }

  template<typename Scalar, int Options>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateCompliance(const Scalar & compliance)
  {
    const Eigen::Index total_constraint_size = constraintDim();
    updateCompliance(Vector::Constant(total_constraint_size, compliance));
  }

  template<typename Scalar, int Options>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::computeDelassusBlock()
  {
    // delassus_block.setZero();
    const auto total_constraint_size = constraintDim();
    const auto UtopRight = U.topRightCorner(total_constraint_size, nv);
    const auto Dtail = D.tail(nv);

    // // Upper left triangular part of U
    //   for (Eigen::Index j = total_constraint_size - 1; j >= 0; --j)
    //   {
    //     const Eigen::Index slice_dim = nv;
    //     typedef Eigen::Map<Vector,EIGEN_DEFAULT_ALIGN_BYTES> MapVector;
    //     MapVector DUt_partial = MapVector(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar,slice_dim,1));

    //     DUt_partial.noalias() =
    //       UtopRight.row(j).transpose().cwiseProduct(Dtail);
    //     for (Eigen::Index _i = j; _i >= 0; _i--)
    //     {
    //       delassus_block(_i, j) = UtopRight.row(_i).dot(DUt_partial);
    //     }
    //   }

    // typedef Eigen::Map<RowMatrix> MapRowMatrix;
    // MapRowMatrix OSIMinv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, total_constraint_size,
    // nv)); OSIMinv.noalias() = UtopRight * Dtail.asDiagonal(); delassus_block.noalias() = OSIMinv
    // * UtopRight.transpose();

    delassus_block.noalias() = (UtopRight * Dtail.asDiagonal()) * UtopRight.transpose();
  }

  template<typename Scalar, int Options>
  template<typename VectorLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateDamping(
    const Eigen::MatrixBase<VectorLike> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
    damping = vec;
    const auto constraint_size = constraintDim();

    auto UTopLeft = U.topLeftCorner(constraint_size, constraint_size);
    UTopLeft.setIdentity();

    // Upper left triangular part of U
    for (Eigen::Index j = constraint_size - 1; j >= 0; --j)
    {
      const Eigen::Index slice_dim = constraint_size - j - 1;

      typedef Eigen::Map<Vector, EIGEN_DEFAULT_ALIGN_BYTES> MapVector;
      MapVector DUt_partial = MapVector(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, slice_dim, 1));
      DUt_partial.noalias() =
        U.row(j).segment(j + 1, slice_dim).transpose().cwiseProduct(D.segment(j + 1, slice_dim));

      D[j] = -delassus_block(j, j) - damping[j] - compliance[j]
             - U.row(j).segment(j + 1, slice_dim).dot(DUt_partial);

      assert(
        check_expression_if_real<Scalar>(D[j] != Scalar(0))
        && "The diagonal element is equal to zero.");
      Dinv[j] = Scalar(1) / D[j];

      for (Eigen::Index _i = j - 1; _i >= 0; _i--)
      {
        U(_i, j) =
          (-delassus_block(_i, j) - U.row(_i).segment(j + 1, slice_dim).dot(DUt_partial)) * Dinv[j];
      }
    }
  }

  template<typename Scalar, int Options>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::updateDamping(const Scalar & mu)
  {
    //      PINOCCHIO_CHECK_INPUT_ARGUMENT(check_expression_if_real<Scalar>(mu >= 0), "mu should be
    //      positive.");
    updateDamping(Vector::Constant(constraintDim(), mu));
  }

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::solveInPlace(
    const Eigen::MatrixBase<MatrixLike> & mat_) const
  {
    auto & mat = mat_.const_cast_derived();

    Uiv(mat);
    mat.array().colwise() *= Dinv.array();
    Utiv(mat);
  }

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::solve(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    Matrix res(mat);
    solveInPlace(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
  ContactCholeskyDecompositionTpl<Scalar, Options>
  ContactCholeskyDecompositionTpl<Scalar, Options>::getMassMatrixChoeslkyDecomposition(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    const DataTpl<S1, O1, JointCollectionTpl> & data) const
  {
    typedef ContactCholeskyDecompositionTpl<Scalar, Options> ReturnType;
    ReturnType res(model, data);

    res.D = D.tail(nv);
    res.Dinv = Dinv.tail(nv);
    res.U = U.bottomRightCorner(nv, nv);

    return res;
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UvAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat_)
      {
        auto & mat = mat_.const_cast_derived();

        assert(mat.rows() == chol.size() && "The input matrix is of wrong size");

        for (Eigen::Index col_id = 0; col_id < mat_.cols(); ++col_id)
          UvAlgo<typename MatrixLike::ColXpr>::run(chol, mat.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UvAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec_)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        auto & vec = vec_.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");
        const Eigen::Index total_constraint_size = chol.constraintDim();

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::Index k = 0; k < total_constraint_size; ++k)
        {
          const Eigen::Index slice_dim = chol.size() - k - 1;
          vec[k] += chol.U.row(k).tail(slice_dim).dot(vec_.tail(slice_dim));
        }

        for (Eigen::Index k = total_constraint_size; k <= chol.size() - 2; ++k)
          vec[k] += chol.U.row(k)
                      .segment(k + 1, chol.nv_subtree_fromRow[k] - 1)
                      .dot(vec.segment(k + 1, chol.nv_subtree_fromRow[k] - 1));
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Uv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UvAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UtvAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = mat.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::Index col_id = 0; col_id < mat_.cols(); ++col_id)
          UtvAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UtvAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = vec.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");
        const Eigen::Index total_constraint_size = chol.constraintDim();

        for (Eigen::Index k = chol.size() - 2; k >= total_constraint_size; --k)
          vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1) +=
            chol.U.row(k).segment(k + 1, chol.nv_subtree_fromRow[k] - 1).transpose() * vec_[k];

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::Index k = total_constraint_size - 1; k >= 0; --k)
        {
          const Eigen::Index slice_dim = chol.size() - k - 1;
          vec_.tail(slice_dim) += chol.U.row(k).tail(slice_dim).transpose() * vec_[k];
        }
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Utv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UtvAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UivAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = mat.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::Index col_id = 0; col_id < mat_.cols(); ++col_id)
          UivAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UivAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = vec.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");

        const Eigen::Index total_constraint_size = chol.constraintDim();
        for (Eigen::Index k = chol.size() - 2; k >= total_constraint_size; --k)
          vec_[k] -= chol.U.row(k)
                       .segment(k + 1, chol.nv_subtree_fromRow[k] - 1)
                       .dot(vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1));

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::Index k = total_constraint_size - 1; k >= 0; --k)
        {
          const Eigen::Index slice_dim = chol.size() - k - 1;
          vec_[k] -= chol.U.row(k).tail(slice_dim).dot(vec_.tail(slice_dim));
        }
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Uiv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UivAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime>
    struct UtivAlgo
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<MatrixLike> & mat)
      {
        MatrixLike & mat_ = mat.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          mat.rows() == chol.size(), "The input matrix is of wrong size");

        for (Eigen::Index col_id = 0; col_id < mat_.cols(); ++col_id)
          UtivAlgo<typename MatrixLike::ColXpr>::run(chol, mat_.col(col_id));
      }
    };

    template<typename VectorLike>
    struct UtivAlgo<VectorLike, 1>
    {
      template<typename Scalar, int Options>
      static void run(
        const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
        const Eigen::MatrixBase<VectorLike> & vec)
      {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike)
        VectorLike & vec_ = vec.const_cast_derived();

        PINOCCHIO_CHECK_INPUT_ARGUMENT(
          vec.size() == chol.size(), "The input vector is of wrong size");
        const Eigen::Index total_constraint_size = chol.constraintDim();

        // TODO: exploit the Sparsity pattern of the first rows of U
        for (Eigen::Index k = 0; k < total_constraint_size; ++k)
        {
          const Eigen::Index slice_dim = chol.size() - k - 1;
          vec_.tail(slice_dim) -= chol.U.row(k).tail(slice_dim).transpose() * vec_[k];
        }

        for (Eigen::Index k = total_constraint_size; k <= chol.size() - 2; ++k)
          vec_.segment(k + 1, chol.nv_subtree_fromRow[k] - 1) -=
            chol.U.row(k).segment(k + 1, chol.nv_subtree_fromRow[k] - 1).transpose() * vec_[k];
      }
    };
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixLike>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::Utiv(
    const Eigen::MatrixBase<MatrixLike> & mat) const
  {
    details::UtivAlgo<MatrixLike>::run(*this, mat.const_cast_derived());
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::matrix() const
  {
    Matrix res(size(), size());
    matrix(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::matrix(
    const Eigen::MatrixBase<MatrixType> & res_) const
  {
    auto & res = res_.const_cast_derived();
    res.noalias() = U * D.asDiagonal() * U.transpose();
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::inverse() const
  {
    Matrix res(size(), size());
    inverse(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1>
  bool ContactCholeskyDecompositionTpl<Scalar, Options>::operator==(
    const ContactCholeskyDecompositionTpl<S1, O1> & other) const
  {
    bool is_same = true;

    if (nv != other.nv)
      return false;

    if (
      D.size() != other.D.size() || Dinv.size() != other.Dinv.size() || U.rows() != other.U.rows()
      || U.cols() != other.U.cols())
      return false;

    is_same &= (D == other.D);
    is_same &= (Dinv == other.Dinv);
    is_same &= (U == other.U);

    is_same &= (parents_fromRow == other.parents_fromRow);
    is_same &= (nv_subtree_fromRow == other.nv_subtree_fromRow);
    //        is_same &= (rowise_sparsity_pattern == other.rowise_sparsity_pattern);

    return is_same;
  }

  template<typename Scalar, int Options>
  template<typename S1, int O1>
  bool ContactCholeskyDecompositionTpl<Scalar, Options>::operator!=(
    const ContactCholeskyDecompositionTpl<S1, O1> & other) const
  {
    return !(*this == other);
  }

  namespace details
  {

    template<typename Scalar, int Options, typename VectorLike>
    PINOCCHIO_DONT_INLINE VectorLike & inverseAlgo(
      const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
      const Eigen::Index col,
      const Eigen::MatrixBase<VectorLike> & vec_)
    {
      EIGEN_STATIC_ASSERT_VECTOR_ONLY(VectorLike);

      auto & vec = vec_.const_cast_derived();
      const Eigen::Index & chol_dim = chol.size();
      PINOCCHIO_CHECK_INPUT_ARGUMENT(col < chol_dim && col >= 0);
      PINOCCHIO_CHECK_INPUT_ARGUMENT(vec.size() == chol_dim);

      const auto & nvt = chol.nv_subtree_fromRow;

      const Eigen::Index last_col =
        std::min(col - 1, chol_dim - 2); // You can start from nv-2 (no child in nv-1)
      vec[col] = Scalar(1);
      vec.tail(chol_dim - col - 1).setZero();

      // TODO: exploit the sparsity pattern of the first rows of U
      for (Eigen::Index k = last_col; k >= 0; --k)
      {
        const Eigen::Index nvt_max = std::min(col - k, nvt[k] - 1);
        const auto U_row = chol.U.row(k);
        vec[k] = -U_row.segment(k + 1, nvt_max).dot(vec.segment(k + 1, nvt_max));
        //  if(k >= chol_constraint_size)
        //  {
        //    vec[k] = -U_row.segment(k+1,nvt_max).dot(vec.segment(k+1,nvt_max));
        //  }
        //  else
        //  {
        //    typedef ContactCholeskyDecompositionTpl<Scalar, Options> ContactCholeskyDecomposition;
        //    typedef typename ContactCholeskyDecomposition::SliceVector SliceVector;
        //    typedef typename ContactCholeskyDecomposition::Slice Slice;
        //    const SliceVector & slice_vector = chol.rowise_sparsity_pattern[(size_t)k];

        //    const Slice & slice_0 = slice_vector[0];
        //    assert(slice_0.first_index == k);
        //    Eigen::Index last_index1 = slice_0.first_index + slice_0.size;
        //    const Eigen::Index last_index2 = k + nvt_max;
        //    Eigen::Index slice_dim = std::min(last_index1,last_index2) - k;
        //    vec[k] =
        //    -U_row.segment(slice_0.first_index+1,slice_dim-1).dot(vec.segment(slice_0.first_index+1,slice_dim-1));

        //    typename SliceVector::const_iterator slice_it = slice_vector.begin()++;
        //    for(;slice_it != slice_vector.end(); ++slice_it)
        //    {
        //      const Slice & slice = *slice_it;
        //      last_index1 = slice.first_index + slice.size;
        //      slice_dim = std::min(last_index1,last_index2+1) - slice.first_index;
        //      if(slice_dim <= 0) break;

        //      vec[k] -=
        //      U_row.segment(slice.first_index,slice_dim).dot(vec_.segment(slice.first_index,slice_dim));
        //    }
        //  }
      }

      vec.head(col + 1).array() *= chol.Dinv.head(col + 1).array();

      for (Eigen::Index k = 0; k < col + 1; ++k) // You can stop one step before nv.
      {
        const Eigen::Index nvt_max = nvt[k] - 1;
        vec.segment(k + 1, nvt_max) -= chol.U.row(k).segment(k + 1, nvt_max).transpose() * vec[k];
      }

      return vec;
    }
  } // namespace details

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::inverse(
    const Eigen::MatrixBase<MatrixType> & res_) const
  {
    auto & res = res_.const_cast_derived();
    PINOCCHIO_CHECK_INPUT_ARGUMENT(res.rows() == size());
    PINOCCHIO_CHECK_INPUT_ARGUMENT(res.cols() == size());

    for (Eigen::Index col_id = 0; col_id < size(); ++col_id)
      details::inverseAlgo(*this, col_id, res.col(col_id));

    res.template triangularView<Eigen::StrictlyLower>() =
      res.transpose().template triangularView<Eigen::StrictlyLower>();
  }

  template<typename Scalar, int Options>
  const typename ContactCholeskyDecompositionTpl<Scalar, Options>::EigenStorageVector::ConstMapType
  ContactCholeskyDecompositionTpl<Scalar, Options>::getCompliance() const
  {
    return compliance_storage.const_map();
  }

  template<typename Scalar, int Options>
  const typename ContactCholeskyDecompositionTpl<Scalar, Options>::EigenStorageVector::ConstMapType
  ContactCholeskyDecompositionTpl<Scalar, Options>::getDamping() const
  {
    return damping_storage.const_map();
  }

  template<typename Scalar, int Options>
  Eigen::Index ContactCholeskyDecompositionTpl<Scalar, Options>::size() const
  {
    return D.size();
  }

  template<typename Scalar, int Options>
  Eigen::Index ContactCholeskyDecompositionTpl<Scalar, Options>::constraintDim() const
  {
    return size() - nv;
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::getInverseOperationalSpaceInertiaMatrix(
    bool enforce_symmetry) const
  {
    Matrix res(constraintDim(), constraintDim());
    getInverseOperationalSpaceInertiaMatrix(res, enforce_symmetry);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::getInverseOperationalSpaceInertiaMatrix(
    const Eigen::MatrixBase<MatrixType> & res, bool enforce_symmetry) const
  {
    const auto U1 = U.topLeftCorner(constraintDim(), constraintDim());

    const auto dim = constraintDim();
    typedef Eigen::Map<RowMatrix, EIGEN_DEFAULT_ALIGN_BYTES> MapRowMatrix;
    MapRowMatrix OSIMinv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, dim, dim));

    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
    MatrixType & res_ = res.const_cast_derived();
    OSIMinv.noalias() = D.head(dim).asDiagonal() * U1.adjoint();
    res_.noalias() = -U1 * OSIMinv;
    if (enforce_symmetry)
      enforceSymmetry(res_);
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::DelassusCholeskyExpression
  ContactCholeskyDecompositionTpl<Scalar, Options>::getDelassusCholeskyExpression() const
  {
    return DelassusCholeskyExpression(*this);
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::getOperationalSpaceInertiaMatrix() const
  {
    Matrix res(constraintDim(), constraintDim());
    getOperationalSpaceInertiaMatrix(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::getOperationalSpaceInertiaMatrix(
    const Eigen::MatrixBase<MatrixType> & res_) const
  {
    auto & res = res_.const_cast_derived();
    //        typedef typename RowMatrix::ConstBlockXpr ConstBlockXpr;
    const auto U1 =
      U.topLeftCorner(constraintDim(), constraintDim()).template triangularView<Eigen::UnitUpper>();

    const auto dim = constraintDim();
    typedef Eigen::Map<RowMatrix, EIGEN_DEFAULT_ALIGN_BYTES> MapRowMatrix;
    MapRowMatrix OSIMinv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, dim, dim));

    typedef Eigen::Map<Matrix, EIGEN_DEFAULT_ALIGN_BYTES> MapMatrix;
    MapMatrix U1inv = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, dim, dim));

    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
    U1inv.setIdentity();
    U1.solveInPlace(U1inv); // TODO: implement Sparse Inverse
    OSIMinv.noalias() = -U1inv.adjoint() * Dinv.head(dim).asDiagonal();
    res.noalias() = OSIMinv * U1inv;
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();
  }

  template<typename Scalar, int Options>
  typename ContactCholeskyDecompositionTpl<Scalar, Options>::Matrix
  ContactCholeskyDecompositionTpl<Scalar, Options>::getInverseMassMatrix() const
  {
    Matrix res(nv, nv);
    getInverseMassMatrix(res);
    return res;
  }

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::getInverseMassMatrix(
    const Eigen::MatrixBase<MatrixType> & res_) const
  {
    auto & res = res_.const_cast_derived();
    //        typedef typename RowMatrix::ConstBlockXpr ConstBlockXpr;
    const auto U4 = U.bottomRightCorner(nv, nv).template triangularView<Eigen::UnitUpper>();

    typedef Eigen::Map<RowMatrix, EIGEN_DEFAULT_ALIGN_BYTES> MapRowMatrix;
    MapRowMatrix Minv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, nv, nv));

    typedef Eigen::Map<Matrix, EIGEN_DEFAULT_ALIGN_BYTES> MapMatrix;
    MapMatrix U4inv = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, nv, nv));

    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
    U4inv.setIdentity();
    U4.solveInPlace(U4inv); // TODO: implement Sparse Inverse
    Minv.noalias() = U4inv.adjoint() * Dinv.tail(nv).asDiagonal();
    res.noalias() = Minv * U4inv;
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();
  }

  template<typename Scalar, int Options>
  template<typename MatrixType>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::getJMinv(
    const Eigen::MatrixBase<MatrixType> & res_) const
  {
    PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
    auto & res = res_.const_cast_derived();
    const auto U4 = U.bottomRightCorner(nv, nv).template triangularView<Eigen::UnitUpper>();
    auto U2 = U.topRightCorner(constraintDim(), nv);

    typedef Eigen::Map<Matrix, EIGEN_DEFAULT_ALIGN_BYTES> MapMatrix;
    MapMatrix U4inv = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, nv, nv));

    U4inv.setIdentity();
    U4.solveInPlace(U4inv); // TODO: implement Sparse Inverse
    res.noalias() = U2 * U4inv;
    PINOCCHIO_EIGEN_MALLOC_ALLOWED();
  }

  template<typename Scalar, int Options>
  template<
    typename S1,
    int O1,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    class ConstraintModelAllocator,
    class ConstraintData,
    class ConstraintDataAllocator>
  void ContactCholeskyDecompositionTpl<Scalar, Options>::compute(
    const ModelTpl<S1, O1, JointCollectionTpl> & model,
    DataTpl<S1, O1, JointCollectionTpl> & data,
    const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
    const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
    const S1 mu)
  {
    compute(
      model, data, constraint_models, constraint_datas, Vector::Constant(constraintDim(), mu));
  }

  template<typename Scalar, int Options>
  std::size_t ContactCholeskyDecompositionTpl<Scalar, Options>::sizeInBytes() const
  {
    return U_storage.sizeInBytes() + D_storage.sizeInBytes() + Dinv_storage.sizeInBytes()
           + compliance_storage.sizeInBytes() + damping_storage.sizeInBytes()
           + delassus_block_storage.sizeInBytes() + pinocchio::sizeInBytes(parents_fromRow)
           + pinocchio::sizeInBytes(nv_subtree_fromRow)
      // + pinocchio::sizeInBytes(rowise_sparsity_pattern)
      ;
  }

  PINOCCHIO_COMPILER_DIAGNOSTIC_POP
} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_constraint_cholesky_hxx__
