//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_algorithm_contact_cholesky_fwd_hpp__
#define __pinocchio_algorithm_contact_cholesky_fwd_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/math/matrix-block.hpp"
#include "pinocchio/math/triangular-matrix.hpp"
#include "pinocchio/container/eigen-storage.hpp"
#include "pinocchio/utils/std-vector.hpp"

#include "pinocchio/algorithm/fwd.hpp"

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

  template<typename _ContactCholeskyDecomposition>
  struct DelassusCholeskyExpressionTpl;

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
  struct PINOCCHIO_DLLAPI
    PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
      ContactCholeskyDecompositionTpl
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef pinocchio::Index Index;
    typedef _Scalar Scalar;
    enum
    {
      LINEAR = 0,
      ANGULAR = 3,
      Options = _Options
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
    typedef typename PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(Matrix) RowMatrix;

    typedef EigenStorageTpl<Vector> EigenStorageVector;
    typedef EigenStorageTpl<Matrix> EigenStorageMatrix;
    typedef EigenStorageTpl<RowMatrix> EigenStorageRowMatrix;

    typedef Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1, Options> EigenIndexVector;
    typedef typename PINOCCHIO_ALIGNED_STD_VECTOR(EigenIndexVector) VectorOfEigenIndexVector;
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
    ContactCholeskyDecompositionTpl();

    ///
    /// \brief Constructor from a model.
    ///
    /// \param[in] model Model of the kinematic tree.
    /// \param[in] data Data associated with the kinematic tree
    ///
    template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
    explicit ContactCholeskyDecompositionTpl(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const DataTpl<S1, O1, JointCollectionTpl> & data);

    ///
    /// \brief Constructor from a model and a collection of RigidConstraintModel objects.
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
      const std::vector<ConstraintData, ConstraintDataAllocator> & constraint_datas);

    ///
    /// \brief Copy constructor
    ///
    /// \param[in] other ContactCholeskyDecompositionTpl to copy
    ///
    ContactCholeskyDecompositionTpl(const ContactCholeskyDecompositionTpl & other);

    ///
    /// \brief Copy operator
    ///
    /// \param[in] other ContactCholeskyDecompositionTpl to copy
    ///
    ContactCholeskyDecompositionTpl & operator=(const ContactCholeskyDecompositionTpl & other);

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
    ///  \brief Internal memory allocation.
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
    void resize(
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
      const S1 mu = S1(0.),
      bool use_explicit_delassus = false);
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
      const Eigen::MatrixBase<VectorLike> & mus,
      bool use_explicit_delassus = false);

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
    const typename EigenStorageVector::MapType getCompliance() const;

    ///
    /// \brief Update the damping terms on the upper left block part of the KKT matrix. The damping
    /// terms should be all positives.
    ///
    /// \param[in] mus Vector of positive regularization factor allowing to enforce the definite
    /// property of the KKT matrix.
    ///
    template<typename VectorLike>
    void
    updateDamping(const Eigen::MatrixBase<VectorLike> & mus, bool use_explicit_delassus = false);

    ///
    /// \brief Update the damping term on the upper left block part of the KKT matrix. The damping
    /// terms should be all positives.
    ///
    /// \param[in] mu Regularization factor allowing to enforce the definite property of the KKT
    /// matrix.
    ///
    void updateDamping(const Scalar & mu, bool use_explicit_delassus = false);

    void computedelassus_blockFromU();

    template<typename VectorLike>
    void updateDampingdelassus_block(const Eigen::MatrixBase<VectorLike> & mus);

    ///
    /// \brief Returns the current damping vector.
    ///
    const typename EigenStorageVector::MapType getDamping() const;

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
    EigenStorageVector D_storage;
    typename EigenStorageVector::RefMapType D;
    EigenStorageVector Dinv_storage;
    typename EigenStorageVector::RefMapType Dinv;
    EigenStorageRowMatrix U_storage;
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
    bool operator==(const ContactCholeskyDecompositionTpl<S1, O1> & other) const;

    template<typename S1, int O1>
    bool operator!=(const ContactCholeskyDecompositionTpl<S1, O1> & other) const;
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

  protected:
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

    EigenStorageRowMatrix delassus_block_storage;
    typename EigenStorageRowMatrix::RefMapType delassus_block;
  };

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_contact_cholesky_fwd_hpp__
