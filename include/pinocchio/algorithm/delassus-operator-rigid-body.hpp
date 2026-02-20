//
// Copyright (c) 2024-2026 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__

#include "pinocchio/math/block-diagonal-matrix.hpp"

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/constraints/utils.hpp"
#include "pinocchio/container/eigen-storage.hpp"

#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/utils/template-template-parameter.hpp"
#include "pinocchio/utils/size-in-bytes.hpp"

namespace pinocchio
{

  template<
    typename _Scalar,
    int _Options,
    template<typename, int> class JointCollectionTpl,
    class _ConstraintModel,
    template<typename T> class Holder>
  struct traits<DelassusOperatorRigidBodySystemsTpl<
    _Scalar,
    _Options,
    JointCollectionTpl,
    _ConstraintModel,
    Holder>>
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    static constexpr int RowsAtCompileTime = Eigen::Dynamic;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> DenseMatrix;
    typedef DenseMatrix Matrix;

    typedef EigenStorageTpl<Vector> EigenStorageVector;
    typedef BlockDiagonalMatrixTpl<Scalar, Options> BlockDiagonalMatrix;

    typedef ModelTpl<Scalar, Options, JointCollectionTpl> Model;
    typedef typename Model::Data Data;

    typedef _ConstraintModel ConstraintModel;
    typedef typename helper::remove_holder<ConstraintModel>::type InnerConstraintModel;
    typedef typename helper::remove_holder<ConstraintModel>::ref_type ConstraintModelReference;
    typedef
      typename std::remove_reference<ConstraintModelReference>::type ConstraintModelReferenceValue;
    static constexpr bool ConstraintModelIsConst =
      std::is_const<ConstraintModelReferenceValue>::value;

    typedef
      typename helper::remove_holder<ConstraintModel>::type::ConstraintData InnerConstraintData;
    typedef typename std::conditional<
      helper::is_type_holder<ConstraintModel>::value,
      typename internal::extract_template_template_parameter<ConstraintModel>::template type<
        typename std::
          conditional<ConstraintModelIsConst, const InnerConstraintData, InnerConstraintData>::
            type>,
      InnerConstraintData>::type ConstraintData;

    typedef std::vector<ConstraintModel> ConstraintModelVector;
    typedef std::vector<ConstraintData> ConstraintDataVector;

    typedef const BlockDiagonalMatrix & getDampingReturnType;
  };

  /// \brief Unsafe version of DelassusOperatorRigidBodySystemsTpl.
  /// Allows to access protected members.
  /// Meant to be used by expert users.
  template<
    typename _Scalar,
    int _Options,
    template<typename, int> class _JointCollectionTpl,
    class _ConstraintModel,
    template<typename T> class _Holder>
  struct Unsafe<DelassusOperatorRigidBodySystemsTpl<
    _Scalar,
    _Options,
    _JointCollectionTpl,
    _ConstraintModel,
    _Holder>>
  {
    typedef DelassusOperatorRigidBodySystemsTpl<
      _Scalar,
      _Options,
      _JointCollectionTpl,
      _ConstraintModel,
      _Holder>
      SafeSelf;
    typedef typename SafeSelf::BlockDiagonalMatrix BlockDiagonalMatrix;

    explicit Unsafe(SafeSelf & self)
    : self(self)
    {
    }

    /// \brief Signal the delassus that updateDecomposition() should be called.
    /// This is typically called after damping or compliance has been updated
    /// so updateSumComplianceDamping is called internally.
    void makeDirty()
    {
      self.updateSumComplianceDamping();
      self.m_solve_in_place_dirty = true;
    }

    /// \brief Getter to the block diagonal damping.
    BlockDiagonalMatrix & damping()
    {
      return self.m_damping;
    }

  protected:
    SafeSelf & self;
  };

  template<
    typename _Scalar,
    int _Options,
    template<typename, int> class _JointCollectionTpl,
    class _ConstraintModel,
    template<typename T> class Holder>
  struct DelassusOperatorRigidBodySystemsTpl
  : DelassusOperatorBase<DelassusOperatorRigidBodySystemsTpl<
      _Scalar,
      _Options,
      _JointCollectionTpl,
      _ConstraintModel,
      Holder>>
  {

    typedef DelassusOperatorRigidBodySystemsTpl Self;
    typedef DelassusOperatorBase<Self> Base;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    typedef typename traits<Self>::Vector Vector;
    typedef typename traits<Self>::DenseMatrix DenseMatrix;
    typedef typename traits<Self>::EigenStorageVector EigenStorageVector;
    typedef typename traits<Self>::BlockDiagonalMatrix BlockDiagonalMatrix;

    typedef typename traits<Self>::Model Model;
    typedef Holder<const Model> ModelHolder;
    typedef typename traits<Self>::Data Data;
    typedef Holder<Data> DataHolder;

    typedef typename Data::Force Force;
    typedef typename Data::VectorXs VectorXs;
    typedef std::vector<Force> ForceVector;

    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::InnerConstraintModel InnerConstraintModel;
    typedef typename traits<Self>::ConstraintModelVector ConstraintModelVector;
    typedef Holder<const ConstraintModelVector> ConstraintModelVectorHolder;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::InnerConstraintData InnerConstraintData;
    typedef typename traits<Self>::ConstraintDataVector ConstraintDataVector;
    typedef Holder<const ConstraintDataVector> ConstraintDataVectorHolder;

    /// \brief Cast this class to its unsafe version.
    Unsafe<Self> unsafe()
    {
      return Unsafe<Self>(*this);
    }
    friend struct Unsafe<Self>;

    /// \brief Default constructor from model, data, constraint_models and constraint_datas.
    DelassusOperatorRigidBodySystemsTpl(
      const ModelHolder & model_ref,
      const DataHolder & data_ref,
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref,
      const Scalar min_damping_value = 0)
    : Base()
    , m_size(residualSize(helper::get_ref(constraint_models_ref)))
    , m_min_damping_value(min_damping_value)
    , m_model_ref(model_ref)
    , m_data_ref(data_ref)
    , m_constraint_models_ref(constraint_models_ref)
    , m_constraint_datas_ref(constraint_datas_ref)
    , m_internal_data(helper::get_ref(model_ref))
    , m_solve_in_place_dirty(true)
    , m_damping(VectorXs::Constant(m_size, min_damping_value).asDiagonal())
    , m_compliance_storage(m_size)
    , m_compliance(m_compliance_storage.map())
    , m_sum_compliance_damping(VectorXs::Constant(m_size, min_damping_value).asDiagonal())
    , m_sum_compliance_damping_inverse(VectorXs::Constant(m_size, min_damping_value).asDiagonal())
    {
      assert(model().check(data()) && "data is not consistent with model.");
      PINOCCHIO_CHECK_ARGUMENT_SIZE(
        constraint_models().size(), constraint_datas().size(),
        "The sizes of contact vector models and contact vector datas are not the same.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        min_damping_value >= Scalar(0) && "The damping value should be positive.");

      rebuild(model_ref, data_ref, constraint_models_ref, constraint_datas_ref);
    }

    /// \brief Update the constraint model and data vectors, and resize the internal quantities.
    ///
    /// \param[in] constraint_models_ref Vector of constraint models
    /// \param[in] constraint_datas_ref Vector of constraint datas
    ///
    void rebuild(
      const ModelHolder & model_ref,
      const DataHolder & data_ref,
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref);

    ///
    /// \brief Update the intermediate computations before calling solveInPlace or operator*
    ///
    /// \param[in] apply_on_the_right If true, this will update the quantities involved in the
    /// applyOnTheRight method
    /// \param[in] solve_in_place If true, this will update the quantities involved in the
    /// solveInPlace method
    ///
    /// \remarks By activating or deactivating apply_on_the_right and solve_in_place, this enables
    /// to lower the quantities updated to the minimum, helping to save time overall.
    /// This method assumes the fields data.oMi, data.lMi and data.J have been computed.
    /// This is typically done by calling `computeJointJacobians` or `aba` in `Convention::WORLD`.
    ///
    void compute(bool apply_on_the_right = true, bool solve_in_place = true)
    {
      compute_or_update_decomposition(apply_on_the_right, solve_in_place);
    }

    /// \brief Stores the product delassus * x in res.
    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x,
      const Eigen::MatrixBase<MatrixOut> & res,
      bool with_damping = true) const;

    /// \brief Update the numerical damping of the delassus from a vector.
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & damping_vector)
    {
      m_damping = damping_vector.asDiagonal();
      updateSumComplianceDamping();
    }

    /// \brief Update the numerical damping of the delassus from a scalar.
    void updateDamping(const Scalar & mu)
    {
      updateDamping(Vector::Constant(size(), mu));
    }

    /// \brief Update numerical damping by copying an input block diagonal matrix.
    template<int OtherOptions, std::size_t OtherAlignment>
    void updateDamping(const BlockDiagonalMatrixTpl<Scalar, OtherOptions, OtherAlignment> &
                         block_diagonal_damping_matrix)
    {
      if (&block_diagonal_damping_matrix == &m_damping)
        return;

      m_damping = block_diagonal_damping_matrix;
      updateSumComplianceDamping();
    }

    /// \brief Update numerical damping by moving an input block diagonal matrix.
    template<int OtherOptions, std::size_t OtherAlignment>
    void updateDamping(
      BlockDiagonalMatrixTpl<Scalar, OtherOptions, OtherAlignment> && block_diagonal_damping_matrix)
    {
      if (&block_diagonal_damping_matrix == &m_damping)
        return;

      m_damping = std::move(block_diagonal_damping_matrix);
      updateSumComplianceDamping();
    }

    /// \brief Update the physical compliance of the delassus from a vector.
    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance_vector)
    {
      m_compliance = compliance_vector;
      updateSumComplianceDamping();
    }

    /// \brief Update the physical compliance of the delassus from a scalar.
    void updateCompliance(const Scalar & compliance_value)
    {
      updateCompliance(Vector::Constant(size(), compliance_value));
    }

    ///
    /// \brief Update the decomposition after a call to updateDamping or updateCompliance.
    ///
    /// \remarks isDirty() allows to retrieve the current status of the decomposition.
    ///
    void updateDecomposition()
    {
      compute(false, true);
    }

    /// \brief Returns true if updateDecomposition() needs to be called in order to call
    /// solveInPlace.
    bool isDirty() const
    {
      return m_solve_in_place_dirty;
    }

    /// \brief Fills input matrix with the dense representation of the Delassus.
    template<typename MatrixType>
    void matrix(
      const Eigen::MatrixBase<MatrixType> & res,
      bool enforce_symmetry = false,
      bool with_damping = true) const
    {
      MatrixType & res_ = res.const_cast_derived();
      typedef Eigen::Map<VectorXs> MapVectorXs;
      MapVectorXs x = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, this->size(), 1));

      for (Eigen::Index i = 0; i < this->size(); ++i)
      {
        x = VectorXs::Unit(this->size(), i);
        this->applyOnTheRight(x, res_.col(i), with_damping);
      }
      if (enforce_symmetry)
      {
        res_ = 0.5 * (res_ + res_.transpose());
      }
    }

    /// \brief Returns the dense representation of the Delassus.
    DenseMatrix matrix(bool enforce_symmetry = false, bool with_damping = true) const
    {
      DenseMatrix res(this->size(), this->size());
      matrix(res, enforce_symmetry, with_damping);
      return res;
    }

    /// \brief Fills input matrix with the dense representation of the Delassus.
    /// The numerical damping is NOT taken into account here.
    template<typename MatrixType>
    void
    undampedMatrix(const Eigen::MatrixBase<MatrixType> & res, bool enforce_symmetry = false) const
    {
      matrix(res, enforce_symmetry, false /*no damping*/);
    }

    /// \brief Returns the dense representation of the Delassus.
    /// The numerical damping is NOT taken into account here.
    DenseMatrix undampedMatrix(bool enforce_symmetry = false) const
    {
      DenseMatrix res(this->size(), this->size());
      matrix(res, enforce_symmetry, false /*no damping*/);
      return res;
    }

    /// \brief solveInPlace operation returning the results of the inverse of the Delassus operator
    /// times the input matrix mat
    ///
    /// \param[in,out] mat Input/output argument containing the right hand side and the result of
    /// the operation
    ///
    /// \warning The parameter is only marked 'const' to make the C++ compiler accept a temporary
    /// expression here. This function will const_cast it, so constness isn't honored here.
    ///
    /// \remarks Even if the method is marked 'const', it will update the internal decomposition if
    /// the Delassus operator is dirty after an update of the damping or compliance values.
    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const;

    /// \brief Returns the current memory footprint of this object in bytes.
    /// \details Sums up the sizes of all internal data members.
    std::size_t sizeInBytes() const
    {
      return m_damping.sizeInBytes() + m_compliance_storage.sizeInBytes()
             + m_sum_compliance_damping.sizeInBytes()
             + m_sum_compliance_damping_inverse.sizeInBytes() + m_internal_data.sizeInBytes();
    }

    /// \brief Returns the number of rows/cols of the Delassus.
    /// The delassus represents a size() x size() linear operator.
    Eigen::Index size() const
    {
      return m_size;
    }

    /// \brief Returns the number of rows of the Delassus.
    Eigen::Index rows() const
    {
      return m_size;
    }

    /// \brief Returns the number of cols of the Delassus.
    Eigen::Index cols() const
    {
      return m_size;
    }

    /// \brief Const getter for model.
    const Model & model() const
    {
      return helper::get_ref(m_model_ref);
    }

    /// \brief Getter for data.
    Data & data()
    {
      return helper::get_ref(m_data_ref);
    }
    ///
    /// \brief Const getter for data.
    const Data & data() const
    {
      return helper::get_ref(m_data_ref);
    }

    /// \brief Const getter of constraint models.
    const ConstraintModelVector & constraint_models() const
    {
      return helper::get_ref(m_constraint_models_ref);
    }

    /// \brief Const getter of constraint datas.
    const ConstraintDataVector & constraint_datas() const
    {
      return helper::get_ref(m_constraint_datas_ref);
    }

    /// \brief Const getter for the numerical damping.
    const BlockDiagonalMatrix & getDamping() const
    {
      return m_damping;
    }

    /// \brief Const getter for the physical compliance.
    typename EigenStorageVector::ConstRefConstMapType getCompliance() const
    {
      return m_compliance_storage.const_map();
    }

    /// \brief Internal data needed for the various passes of the delassus operator.
    struct InternalData
    {
      typedef typename Data::Motion Motion;
      typedef typename Data::Force Force;

      typedef typename std::vector<Motion> MotionVector;
      typedef typename std::vector<Force> ForceVector;

      InternalData(const Model & model)
      : a(std::size_t(model.njoints), Motion::Zero())
      , oa_augmented(std::size_t(model.njoints), Motion::Zero())
      , u_storage(model.nv)
      , u(u_storage.map())
      , ddq_storage(model.nv)
      , ddq(ddq_storage.map())
      , f(std::size_t(model.njoints))
      , of_augmented(std::size_t(model.njoints))
      {
      }

      /// \brief Rebuild from a given model.
      void rebuild(const Model & model)
      {
        a.resize(std::size_t(model.njoints));
        oa_augmented.resize(std::size_t(model.njoints));
        u_storage.resize(model.nv);
        ddq_storage.resize(model.nv);
        f.resize(std::size_t(model.njoints));
        of_augmented.resize(std::size_t(model.njoints));

        assert(u.size() == model.nv);
        assert(ddq.size() == model.nv);
      }

      MotionVector a, oa_augmented;
      EigenStorageVector u_storage;
      typename EigenStorageVector::RefMapType u;
      EigenStorageVector ddq_storage;
      typename EigenStorageVector::RefMapType ddq;
      ForceVector f, of_augmented;

      /// \brief Returns the current memory footprint of this object in bytes.
      /// \details Sums up the sizes of all internal data members.
      std::size_t sizeInBytes() const
      {
        return pinocchio::sizeInBytes(a) + pinocchio::sizeInBytes(oa_augmented)
               + pinocchio::sizeInBytes(u) + pinocchio::sizeInBytes(ddq) + pinocchio::sizeInBytes(f)
               + pinocchio::sizeInBytes(of_augmented);
      }
    };

    /// \brief Const getter for internal data.
    const InternalData & getInternalData() const
    {
      return m_internal_data;
    }

    /// \brief Getter for internal data.
    InternalData & getInternalData()
    {
      return m_internal_data;
    }

    /// \brief AugmentedMassMatrixOperator needed for solveInPlace.
    struct AugmentedMassMatrixOperator
    {
      AugmentedMassMatrixOperator(const DelassusOperatorRigidBodySystemsTpl & delassus_operator)
      : m_self(delassus_operator)
      {
      }

      template<typename MatrixLike>
      void solveInPlace(
        const Eigen::MatrixBase<MatrixLike> & mat, bool reset_joint_force_vector = true) const;

    protected:
      const DelassusOperatorRigidBodySystemsTpl & m_self;
    };

    /// \brief Getter for the AugmentedMassMatrixOperator.
    AugmentedMassMatrixOperator getAugmentedMassMatrixOperator() const
    {
      return AugmentedMassMatrixOperator(*this);
    }

  protected:
    /// \brief Update the sum compliance + damping
    void updateSumComplianceDamping()
    {
      m_sum_compliance_damping = m_damping + m_compliance.asDiagonal();
      m_solve_in_place_dirty = true;
    }

    /// \brief Update the intermediate computations before calling solveInPlace or operator*.
    void compute_or_update_decomposition(bool apply_on_the_right, bool solve_in_place);

    // Holders
    Eigen::Index m_size;
    const Scalar m_min_damping_value;
    ModelHolder m_model_ref;
    DataHolder m_data_ref;
    ConstraintModelVectorHolder m_constraint_models_ref;
    ConstraintDataVectorHolder m_constraint_datas_ref;

    mutable InternalData m_internal_data;
    mutable bool m_solve_in_place_dirty;

    BlockDiagonalMatrix m_damping;
    EigenStorageVector m_compliance_storage;
    typename EigenStorageVector::RefMapType m_compliance;
    BlockDiagonalMatrix m_sum_compliance_damping;
    BlockDiagonalMatrix m_sum_compliance_damping_inverse;
  };

} // namespace pinocchio

#include "pinocchio/algorithm/delassus-operator-rigid-body.hxx"

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
