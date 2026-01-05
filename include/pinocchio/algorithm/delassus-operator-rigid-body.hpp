//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
#define __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__

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
    typename Scalar,
    int Options,
    template<typename, int> class JointCollectionTpl,
    class ConstraintModel,
    template<typename T> class Holder = std::reference_wrapper>
  struct DelassusOperatorRigidBodySystemsTpl;

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

    enum
    {
      Options = _Options,
      RowsAtCompileTime = Eigen::Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> DenseMatrix;
    typedef DenseMatrix Matrix;

    typedef EigenStorageTpl<Vector> EigenStorageVector;

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

    typedef PINOCCHIO_ALIGNED_STD_VECTOR(ConstraintModel) ConstraintModelVector;
    typedef PINOCCHIO_ALIGNED_STD_VECTOR(ConstraintData) ConstraintDataVector;

    typedef const Vector & getDampingReturnType;
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
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef DelassusOperatorRigidBodySystemsTpl Self;
    typedef DelassusOperatorBase<Self> Base;

    typedef typename traits<Self>::Scalar Scalar;
    enum
    {
      Options = traits<Self>::Options
    };

    typedef typename traits<Self>::Vector Vector;
    typedef typename traits<Self>::DenseMatrix DenseMatrix;
    typedef typename traits<Self>::EigenStorageVector EigenStorageVector;

    typedef typename traits<Self>::Model Model;
    typedef Holder<const Model> ModelHolder;
    typedef typename traits<Self>::Data Data;
    typedef Holder<Data> DataHolder;

    typedef typename Data::Force Force;
    typedef typename Data::VectorXs VectorXs;
    typedef PINOCCHIO_ALIGNED_STD_VECTOR(Force) ForceVector;

    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::InnerConstraintModel InnerConstraintModel;
    typedef typename traits<Self>::ConstraintModelVector ConstraintModelVector;
    typedef Holder<const ConstraintModelVector> ConstraintModelVectorHolder;

    typedef typename traits<Self>::ConstraintData ConstraintData;
    typedef typename traits<Self>::InnerConstraintData InnerConstraintData;
    typedef typename traits<Self>::ConstraintDataVector ConstraintDataVector;
    typedef Holder<const ConstraintDataVector> ConstraintDataVectorHolder;

    DelassusOperatorRigidBodySystemsTpl(
      const ModelHolder & model_ref,
      const DataHolder & data_ref,
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref,
      const Scalar min_damping_value = Scalar(1e-8))
    : Base()
    , m_size(
        residualSize(helper::get_ref(constraint_models_ref), helper::get_ref(constraint_datas_ref)))
    , m_min_damping_value(min_damping_value)
    , m_model_ref(model_ref)
    , m_data_ref(data_ref)
    , m_constraint_models_ref(constraint_models_ref)
    , m_constraint_datas_ref(constraint_datas_ref)
    , m_internal_data(helper::get_ref(model_ref))
    , m_solve_in_place_dirty(true)
    , m_damping_storage(m_size)
    , m_damping(m_damping_storage.map())
    , m_compliance_storage(m_size)
    , m_compliance(m_compliance_storage.map())
    , m_sum_compliance_damping_storage(m_size)
    , m_sum_compliance_damping(m_sum_compliance_damping_storage.map())
    , m_sum_compliance_damping_inverse_storage(m_size)
    , m_sum_compliance_damping_inverse(m_sum_compliance_damping_inverse_storage.map())
    {
      assert(model().check(data()) && "data is not consistent with model.");
      PINOCCHIO_CHECK_ARGUMENT_SIZE(
        constraint_models().size(), constraint_datas().size(),
        "The sizes of contact vector models and contact vector datas are not the same.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        min_damping_value >= Scalar(0) && "The damping value should be positive.");

      m_damping.fill(m_min_damping_value);
      m_compliance.setZero();
      update(constraint_models_ref, constraint_datas_ref);
    }

    /// \brief Update the constraint model and data vectors, and resize the internal quantities.
    ///
    /// \param[in] constraint_models_ref Vector of constraint models
    /// \param[in] constraint_datas_ref Vector of constraint datas
    ///
    void update(
      const ConstraintModelVectorHolder & constraint_models_ref,
      const ConstraintDataVectorHolder & constraint_datas_ref);

    template<typename MatrixType>
    void matrix(const Eigen::MatrixBase<MatrixType> & res, bool enforce_symmetry = false) const
    {
      MatrixType & res_ = res.const_cast_derived();
      typedef Eigen::Map<VectorXs> MapVectorXs;
      MapVectorXs x = MapVectorXs(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, this->size(), 1));

      for (Eigen::Index i = 0; i < this->size(); ++i)
      {
        x = VectorXs::Unit(this->size(), i);
        this->applyOnTheRight(x, res_.col(i));
      }
      if (enforce_symmetry)
      {
        res_ = 0.5 * (res_ + res_.transpose());
      }
    }

    DenseMatrix matrix(bool enforce_symmetry = false) const
    {
      DenseMatrix res(this->size(), this->size());
      matrix(res, enforce_symmetry);
      return res;
    }

  protected:
    void compute_or_update_decomposition(bool apply_on_the_right, bool solve_in_place);

  public:
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

    ///
    /// \brief Update the decomposition after a call to updateDamping or updateCompliance.
    ///
    /// \remarks isDirty() allows to retrieve the current status of the decomposition.
    ///
    void updateDecomposition()
    {
      compute(false, true);
    }

  public:
    const Model & model() const
    {
      return helper::get_ref(m_model_ref);
    }

    Data & data()
    {
      return helper::get_ref(m_data_ref);
    }
    const Data & data() const
    {
      return helper::get_ref(m_data_ref);
    }

    const ConstraintModelVector & constraint_models() const
    {
      return helper::get_ref(m_constraint_models_ref);
    }

    const ConstraintDataVector & constraint_datas() const
    {
      return helper::get_ref(m_constraint_datas_ref);
    }

    Eigen::Index size() const
    {
      return m_size;
    }
    Eigen::Index rows() const
    {
      return m_size;
    }
    Eigen::Index cols() const
    {
      return m_size;
    }

    bool isDirty() const
    {
      return m_solve_in_place_dirty;
    }

    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const;

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & damping_vector)
    {
      m_damping = damping_vector;
      updateSumComplianceDamping();
    }

    void updateDamping(const Scalar & mu)
    {
      updateDamping(Vector::Constant(size(), mu));
    }

    typename EigenStorageVector::ConstRefConstMapType getDamping() const
    {
      return m_damping_storage.const_map();
    }

    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance_vector)
    {
      m_compliance = compliance_vector;
      updateSumComplianceDamping();
    }

    void updateCompliance(const Scalar & compliance_value)
    {
      updateCompliance(Vector::Constant(size(), compliance_value));
    }

    typename EigenStorageVector::ConstRefConstMapType getCompliance() const
    {
      return m_compliance_storage.const_map();
    }

    template<typename MatrixLike>
    void updateBarrierHessian(const std::vector<MatrixLike> & blocks)
    {
      PINOCCHIO_UNUSED_VARIABLE(blocks);
      PINOCCHIO_THROW(
        std::runtime_error,
        "updateBarrierHessian not implemented for DelassusOperatorRigidBodySystemsTpl.");
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

    struct InternalData
    {
      typedef typename Data::Motion Motion;
      typedef typename Data::Force Force;

      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(Motion) MotionVector;
      typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(Force) ForceVector;

      InternalData(const Model & model)
      : a(size_t(model.njoints), Motion::Zero())
      , oa_augmented(size_t(model.njoints), Motion::Zero())
      , u(model.nv)
      , ddq(model.nv)
      , f(size_t(model.njoints))
      , of_augmented(size_t(model.njoints))
      {
      }

      MotionVector a, oa_augmented;
      VectorXs u, ddq;
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

    const InternalData & getInternalData() const
    {
      return m_internal_data;
    }

    InternalData & getInternalData()
    {
      return m_internal_data;
    }

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

    AugmentedMassMatrixOperator getAugmentedMassMatrixOperator() const
    {
      return AugmentedMassMatrixOperator(*this);
    }

    /// \brief Returns the current memory footprint of this object in bytes.
    /// \details Sums up the sizes of all internal data members.
    std::size_t sizeInBytes() const
    {
      return m_damping_storage.sizeInBytes() + m_compliance_storage.sizeInBytes()
             + m_sum_compliance_damping_storage.sizeInBytes()
             + m_sum_compliance_damping_inverse_storage.sizeInBytes()
             + m_internal_data.sizeInBytes();
    }

  protected:
    DelassusOperatorRigidBodySystemsTpl & self_const_cast() const
    {
      return const_cast<DelassusOperatorRigidBodySystemsTpl &>(*this);
    }

    void updateSumComplianceDamping() const
    {
      m_sum_compliance_damping = m_damping + m_compliance;
      m_sum_compliance_damping_inverse = m_sum_compliance_damping.cwiseInverse();
      m_solve_in_place_dirty = true;
    }

    // Holders
    Eigen::Index m_size;
    const Scalar m_min_damping_value;
    ModelHolder m_model_ref;
    DataHolder m_data_ref;
    ConstraintModelVectorHolder m_constraint_models_ref;
    ConstraintDataVectorHolder m_constraint_datas_ref;

    mutable InternalData m_internal_data;
    mutable bool m_solve_in_place_dirty;

    EigenStorageVector m_damping_storage;
    typename EigenStorageVector::RefMapType m_damping;
    EigenStorageVector m_compliance_storage;
    typename EigenStorageVector::RefMapType m_compliance;
    EigenStorageVector m_sum_compliance_damping_storage;
    typename EigenStorageVector::RefMapType m_sum_compliance_damping;
    EigenStorageVector m_sum_compliance_damping_inverse_storage;
    typename EigenStorageVector::RefMapType m_sum_compliance_damping_inverse;
  };

} // namespace pinocchio

#include "pinocchio/algorithm/delassus-operator-rigid-body.hxx"

#endif // ifndef __pinocchio_algorithm_delassus_operator_linear_complexity_hpp__
