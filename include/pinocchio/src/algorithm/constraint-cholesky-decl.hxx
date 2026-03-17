//
// Copyright (c) 2019-2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/algorithm/constraint-cholesky.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/algorithm/constraint-cholesky-decl.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  // Forward declaration of algo
  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UtvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UivAlgo;

    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UtivAlgo;

    template<typename Scalar, int Options, typename VectorLike>
    VectorLike & inverseAlgo(
      const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
      const Eigen::Index col,
      const Eigen::MatrixBase<VectorLike> & vec);
  } // namespace details

  ///
  ///  \brief Contact Cholesky decomposition structure. This structure allows
  ///        to compute in a efficient and parsimonious way the Cholesky decomposition
  ///        of the KKT matrix related to the contact dynamics.
  ///        Such a decomposition is usefull when computing both the forward dynamics in contact
  ///        or the related analytical derivatives.
  ///
  ///
  /// \tparam _Scalar Scalar type.
  ///  \tparam _Options Alignment Options of the Eigen objects contained in the data structure.
  ///
  template<typename _Scalar, int _Options>
  struct ContactCholeskyDecompositionTpl
  {

    typedef pinocchio::Index Index;
    typedef _Scalar Scalar;
    static constexpr int LINEAR = 0;
    static constexpr int ANGULAR = 3;
    static constexpr int Options = _Options;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
    typedef typename PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(Matrix) RowMatrix;

    typedef EigenStorageTpl<Vector> EigenStorageVector;
    typedef EigenStorageTpl<Matrix> EigenStorageMatrix;
    typedef EigenStorageTpl<RowMatrix> EigenStorageRowMatrix;

    typedef Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1, Options> EigenIndexVector;
    typedef typename std::vector<EigenIndexVector> VectorOfEigenIndexVector;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector;

    ///@{
    /// \brief Data information related to the Sparsity structure of the Cholesky decompostion
    struct Slice
    {
      Slice(const Eigen::Index & first_index, const Eigen::Index & size)
      : first_index(first_index)
      , size(size)
      {
      }

      Eigen::Index first_index;
      Eigen::Index size;
    };

    typedef DelassusCholeskyExpressionTpl<ContactCholeskyDecompositionTpl>
      DelassusCholeskyExpression;
    friend struct DelassusCholeskyExpressionTpl<ContactCholeskyDecompositionTpl>;

    typedef std::vector<Slice> SliceVector;
    typedef std::vector<SliceVector> VectorOfSliceVector;
    ///@}

    ///
    /// \brief Default constructor
    ///
    explicit ContactCholeskyDecompositionTpl(const Scalar min_damping_value = Scalar(0))
    : D(D_storage.map())
    , Dinv(Dinv_storage.map())
    , U(U_storage.map())
    , compliance(compliance_storage.map())
    , damping(damping_storage.map())
    , sum_compliance_damping(sum_compliance_damping_storage.map())
    , delassus_block(delassus_block_storage.map())
    , decomposition_dirty(true)
    , min_damping_value(min_damping_value)
    {
    }

    ///
    /// \brief Constructor from a model.
    ///
    /// \param[in] model Model of the kinematic tree.
    /// \param[in] data Data associated with the kinematic tree
    ///
    template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
    explicit ContactCholeskyDecompositionTpl(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const DataTpl<S1, O1, JointCollectionTpl> & data,
      const Scalar min_damping_value = Scalar(0));

    ///
    /// \brief Constructor from a model and a collection of ConstraintModel objects.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] data Data associated with the kinematic tree
    /// \param[in] constraint_models Vector of ConstraintModels
    /// \param[in] constraint_datas Vector of constraint datas
    /// information
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    ContactCholeskyDecompositionTpl(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const Scalar min_damping_value = Scalar(0));

    ///
    /// \brief Copy constructor
    ///
    /// \param[in] other ContactCholeskyDecompositionTpl to copy
    ///
    ContactCholeskyDecompositionTpl(const ContactCholeskyDecompositionTpl & other)
    : ContactCholeskyDecompositionTpl(other.min_damping_value)
    {
      *this = other;
    }

    ///
    /// \brief Copy operator
    ///
    /// \param[in] other ContactCholeskyDecompositionTpl to copy
    ///
    ContactCholeskyDecompositionTpl & operator=(const ContactCholeskyDecompositionTpl & other)

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
      sum_compliance_damping_storage = other.sum_compliance_damping_storage;
      delassus_block_storage = other.delassus_block_storage;

      decomposition_dirty = other.decomposition_dirty;
      min_damping_value = other.min_damping_value;

      return *this;
    }

    ///
    ///  \brief Internal memory allocation.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] constraint_models Vector of ConstraintModel
    /// \param[in] constraint_datas Vector of ConstraintData
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    PINOCCHIO_DEPRECATED void allocate(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

    ///
    ///  \brief Rebuild internal memory from a given model and constraints.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] constraint_models Vector of constraint models
    /// \param[in] constraint_datas Vector of constraint datas
    ///
    /// \note This method assumes that the constrained datas are up-to-date.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    void rebuild(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

    ///
    /// \brief Returns the Inverse of the Operational Space Inertia Matrix resulting from the
    /// decomposition.
    ///
    Matrix getInverseOperationalSpaceInertiaMatrix(bool enforce_symmetry = false) const;

    template<typename MatrixType>
    void getInverseOperationalSpaceInertiaMatrix(
      const Eigen::MatrixBase<MatrixType> & res, bool enforce_symmetry = false) const;

    /// \brief Returns the Cholesky decomposition expression associated to the underlying
    /// delassus_block matrix.
    DelassusCholeskyExpression getDelassusCholeskyExpression() const;

    ///
    /// \brief Returns the Operational Space Inertia Matrix resulting from the decomposition.
    ///
    Matrix getOperationalSpaceInertiaMatrix() const;

    template<typename MatrixType>
    void getOperationalSpaceInertiaMatrix(const Eigen::MatrixBase<MatrixType> & res_) const;

    Matrix getInverseMassMatrix() const;

    template<typename MatrixType>
    void getInverseMassMatrix(const Eigen::MatrixBase<MatrixType> & res_) const;

    template<typename MatrixType>
    void getJMinv(const Eigen::MatrixBase<MatrixType> & res_) const;

    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of onstraintModel named constraint_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] constraint_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in] constraint_datas Vector containing the contact data related to the
    /// constraint_models.
    /// \param[in] apply_on_the_right compute quantities related to applyOnTheRight.
    /// \param[in] solve_in_place compute quantities related to solveInPlace.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// This method assumes that the constrained datas are up-to-date.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      bool apply_on_the_right = true,
      bool solve_in_place = true);

    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of ConstraintModel named constraint_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] constraint_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in,out] constraint_datas Vector containing the contact data related to the
    /// constraint_models.
    /// \param[in] mu Positive regularization factor allowing to enforce the definite property of
    /// the KKT matrix.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// This method assumes that the constrained datas are up-to-date.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const S1 mu);
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of onstraintModel named constraint_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] constraint_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in,out] constraint_datas Vector containing the contact data related to the
    /// constraint_models.
    /// \param[in] mu Positive regularization factor allowing to enforce the definite property of
    /// the KKT matrix.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// This method assumes that the constrained datas are up-to-date.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator,
      typename VectorLike>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & constraint_models,
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas,
      const Eigen::MatrixBase<VectorLike> & mus);

    ///
    /// \brief Computes the Cholesky decompostion of the Delassus part. This function should be
    /// called after a call to updateDamping() or updateCompliance().
    ///
    void computeDelassusCholeskyDecomposition();

    ///
    /// \brief Update the compliance terms on the upper left block part of the KKT matrix. The
    /// compliance terms should be all positives.
    ///
    /// \param[in] compliance Vector of physical compliance for the constraints.
    ///
    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance);

    ///
    /// \brief Update the compliance term on the upper left block part of the KKT matrix. The
    /// compliance terms should be all positives.
    ///
    /// \param[in] compliance The physical compliance for the constraints.
    ///
    void updateCompliance(const Scalar & compliance);

    ///
    /// \brief Returns the current compliance vector.
    ///
    const typename EigenStorageVector::ConstMapType getCompliance() const;

    ///
    /// \brief Update the damping terms on the upper left block part of the KKT matrix. The damping
    /// terms should be all positives.
    ///
    /// \param[in] mus Vector of positive regularization factor allowing to enforce the definite
    /// property of the KKT matrix.
    ///
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & mus);

    ///
    /// \brief Update the damping term on the upper left block part of the KKT matrix. The damping
    /// terms should be all positives.
    ///
    /// \param[in] mu Regularization factor allowing to enforce the definite property of the KKT
    /// matrix.
    ///
    void updateDamping(const Scalar & mu);

    ///
    /// \brief Returns the current damping vector.
    ///
    const typename EigenStorageVector::ConstMapType getDamping() const;

    /// \brief Size of the decomposition.
    Eigen::Index size() const;

    /// \brief Returns the total dimension of the constraints contained in the Cholesky
    /// factorization.
    Eigen::Index constraintDim() const;

    ///
    ///  \brief Computes the solution of \f$ A x = b \f$ where *this is the Cholesky decomposition
    /// of A.         "in-place" version of ContactCholeskyDecompositionTpl::solve(b) where the
    /// result is written in b.
    ///        This functions takes as input the vector b, and returns the solution \f$ x = A^-1 b
    ///        \f$.
    ///
    /// \param[inout] mat The right-and-side term which also contains the solution of the linear
    /// system.
    ///
    /// \sa ContactCholeskyDecompositionTpl::solve
    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const;

    ///
    ///  \brief Computes the solution of \f$ A x = b \f$ where *this is the Cholesky decomposition
    /// of A.
    ///        This functions takes as input the vector b, and returns the solution \f$ x = A^-1 b
    ///        \f$.
    ///
    /// \param[inout] mat The right-and-side term.
    ///
    /// \sa ContactCholeskyDecompositionTpl::solveInPlace
    template<typename MatrixLike>
    Matrix solve(const Eigen::MatrixBase<MatrixLike> & mat) const;

    ///
    ///  \brief Retrieves the Cholesky decomposition of the Mass Matrix contained in *this.
    ///
    /// \param[in] model Model of the dynamical system.
    ///
    template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
    ContactCholeskyDecompositionTpl getMassMatrixChoeslkyDecomposition(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const DataTpl<S1, O1, JointCollectionTpl> & data) const;

    ///@{
    /// \brief Vectorwize operations
    template<typename MatrixLike>
    void Uv(const Eigen::MatrixBase<MatrixLike> & mat) const;

    template<typename MatrixLike>
    void Utv(const Eigen::MatrixBase<MatrixLike> & mat) const;

    template<typename MatrixLike>
    void Uiv(const Eigen::MatrixBase<MatrixLike> & mat) const;

    template<typename MatrixLike>
    void Utiv(const Eigen::MatrixBase<MatrixLike> & mat) const;
    ///@}

    /// \brief Returns the matrix resulting from the decomposition
    Matrix matrix() const;

    /// \brief Fill the input matrix with the matrix resulting from the decomposition
    template<typename MatrixType>
    void matrix(const Eigen::MatrixBase<MatrixType> & res) const;

    /// \brief Returns the inverse matrix resulting from the decomposition
    Matrix inverse() const;

    /// \brief Fill the input matrix with the inverse matrix resulting from the decomposition
    template<typename MatrixType>
    void inverse(const Eigen::MatrixBase<MatrixType> & res) const;

    // data
  protected:
    // temporary containing the results of D * U^t
    Vector DUt_storage;
    EigenStorageVector D_storage;
    EigenStorageVector Dinv_storage;
    EigenStorageRowMatrix U_storage;

    /// \brief Inverse of the top left block of U
    mutable Matrix U1inv_storage;
    /// \brief Inverse of the bottom right block of U
    mutable Matrix U4inv_storage;
    mutable RowMatrix OSIMinv_storage, Minv_storage;

  public:
    typename EigenStorageVector::RefMapType D;
    typename EigenStorageVector::RefMapType Dinv;
    typename EigenStorageRowMatrix::RefMapType U;

    ///@{
    /// \brief Friend algorithms
    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UtvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UivAlgo;

    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UtivAlgo;

    // TODO Remove when API is stabilized
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    template<typename Scalar, int Options, typename VectorLike>
    friend VectorLike & details::inverseAlgo(
      const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
      const Eigen::Index col,
      const Eigen::MatrixBase<VectorLike> & vec);
    ///@}

    template<typename S1, int O1>
    bool operator==(const ContactCholeskyDecompositionTpl<S1, O1> & other) const
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

      is_same &= (compliance_storage == other.compliance_storage);
      is_same &= (damping_storage == other.damping_storage);
      is_same &= (sum_compliance_damping_storage == other.sum_compliance_damping_storage);
      is_same &= (delassus_block_storage == other.delassus_block_storage);
      is_same &= (decomposition_dirty == other.decomposition_dirty);
      return is_same;
    }

    template<typename S1, int O1>
    bool operator!=(const ContactCholeskyDecompositionTpl<S1, O1> & other) const
    {
      return !(*this == other);
    }

    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    /// \brief Returns the current memory footprint of this object in bytes.
    /// \details Sums up the sizes of all internal data members.
    std::size_t sizeInBytes() const;

    /// \returns True if the decomposition is dirty, related to change of the damping or compliance
    /// terms.
    /// \remarks compute() will solve the issue by performing a new Cholesky decomposition of the
    /// Delassus block.
    bool isDirty() const
    {
      return decomposition_dirty;
    }

  protected:
    void computeDelassusMatrix();

    void updateSumComplianceDamping();

    EigenIndexVector parents_fromRow;
    EigenIndexVector nv_subtree_fromRow;

    /// \brief Dimension of the tangent of the configuration space of the model
    Eigen::Index nv;

    VectorOfSliceVector rowise_sparsity_pattern;

    /// \brief Store the current value of the physical compliance
    EigenStorageVector compliance_storage;
    typename EigenStorageVector::RefMapType compliance;

    /// \brief Store the current damping value
    EigenStorageVector damping_storage;
    typename EigenStorageVector::RefMapType damping;

    /// \brief Store the sum of compliance and damping vectors
    EigenStorageVector sum_compliance_damping_storage;
    typename EigenStorageVector::RefMapType sum_compliance_damping;

    /// \brief Store the Delassus block
    EigenStorageRowMatrix delassus_block_storage;
    typename EigenStorageRowMatrix::RefMapType delassus_block;

    /// \brief Check if the decomposition is dirty
    bool decomposition_dirty;

    /// \brief Minimum damping value that uninitialized damping is set at.
    Scalar min_damping_value;
  };

} // namespace pinocchio

#ifdef PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #ifndef PINOCCHIO_SKIP_ALGORITHM_CONTACT_CHOLESKY

// Because of a GCC bug we should NEVER define a function that use ContactCholeskyDecompositionTpl
// before doing the explicit template instantiation.
// If we don't take care, GCC will not accept any visibility attribute when declaring the
// explicit template instantiation of the ContactCholeskyDecompositionTpl class.
// The warning message will look like this: type attributes ignored after type is already defined
// [-Wattributes] A minimal code example is added on the PR
// (https://github.com/stack-of-tasks/pinocchio/pull/2469)
namespace pinocchio
{
  extern template struct PINOCCHIO_EXPLICIT_INSTANTIATION_DECLARATION_DLLAPI
    ContactCholeskyDecompositionTpl<context::Scalar, context::Options>;
}

  #endif // PINOCCHIO_SKIP_ALGORITHM_CONTACT_CHOLESKY
#endif   // ifdef PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
