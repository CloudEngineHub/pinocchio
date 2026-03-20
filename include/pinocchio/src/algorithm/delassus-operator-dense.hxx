//
// Copyright (c) 2024-2026 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/algorithm/delassus-operator.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/algorithm/delassus-operator.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  template<
    typename _Scalar,
    int _Options,
    template<typename, auto...> class CholeskyDecompositionTpl>
  struct traits<DelassusOperatorDenseTpl<_Scalar, _Options, CholeskyDecompositionTpl>>
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    static constexpr int RowsAtCompileTime = Eigen::Dynamic;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;

    typedef BlockDiagonalMatrixTpl<Scalar, Options> BlockDiagonalMatrix;

    typedef const BlockDiagonalMatrix & getDampingReturnType;
  };

  /// \brief Unsafe version of DelassusOperatorDenseTpl.
  /// Allows to access protected members.
  /// Meant to be used by expert users.
  template<
    typename _Scalar,
    int _Options,
    template<typename, auto...> class CholeskyDecompositionTpl>
  struct Unsafe<DelassusOperatorDenseTpl<_Scalar, _Options, CholeskyDecompositionTpl>>
  {
    typedef DelassusOperatorDenseTpl<_Scalar, _Options, CholeskyDecompositionTpl> SafeSelf;
    typedef typename SafeSelf::BlockDiagonalMatrix BlockDiagonalMatrix;

    explicit Unsafe(SafeSelf & self)
    : self(self)
    {
    }

    /// \brief Signal the delassus that updateDecomposition() should be called.
    void makeDirty()
    {
      self.m_cholesky_decomposition_dirty = true;
    }

    /// \brief Getter to the block diagonal damping.
    BlockDiagonalMatrix & damping()
    {
      return self.m_damping;
    }

  protected:
    SafeSelf & self;
  };

  /// \brief Operator for a delassus' dense representation.
  template<
    typename _Scalar,
    int _Options,
    template<typename, auto...> class CholeskyDecompositionTpl>
  struct DelassusOperatorDenseTpl
  : DelassusOperatorBase<DelassusOperatorDenseTpl<_Scalar, _Options>>
  {

    typedef _Scalar Scalar;
    typedef DelassusOperatorDenseTpl Self;
    static constexpr int Options = _Options;
    static constexpr int RowsAtCompileTime = traits<DelassusOperatorDenseTpl>::RowsAtCompileTime;

    typedef typename traits<Self>::Matrix Matrix;
    typedef typename traits<Self>::Vector Vector;
    typedef EigenStorageTpl<Matrix> MatrixStorage;
    typedef EigenStorageTpl<Vector> VectorStorage;
    typedef typename MatrixStorage::RefMapType MatrixStorageRefMapType;
    typedef typename VectorStorage::RefMapType VectorStorageRefMapType;
    typedef typename traits<Self>::BlockDiagonalMatrix BlockDiagonalMatrix;
    typedef CholeskyDecompositionTpl<Eigen::Ref<Matrix>> CholeskyDecomposition;
    typedef DelassusOperatorBase<Self> Base;

    /// \brief Cast this class to its unsafe version.
    Unsafe<Self> unsafe()
    {
      return Unsafe<Self>(*this);
    }
    friend struct Unsafe<Self>;

    /// \brief Default constructor.
    DelassusOperatorDenseTpl()
    : Base()
    , m_cholesky_decomposition(m_cholesky_decomposition_data)
    , m_cholesky_decomposition_dirty(true)
    {
    }

    /// \brief Constructor from a given matrix.
    /// \note The constructor does not compute the cholesky decomposition of the delassus.
    template<typename MatrixDerived>
    explicit DelassusOperatorDenseTpl(const Eigen::MatrixBase<MatrixDerived> & mat)
    : DelassusOperatorDenseTpl()
    {
      rebuild(mat);
    }

    /// \brief Constructor from a DelassusCholeskyExpressionTpl.
    /// \note The constructor does not compute the cholesky decomposition of the delassus.
    template<typename ConstraintCholeskyDecomposition>
    explicit DelassusOperatorDenseTpl(
      const DelassusCholeskyExpressionTpl<ConstraintCholeskyDecomposition> & delassus_expression,
      const bool enforce_symmetry = false)
    : DelassusOperatorDenseTpl()
    {
      rebuild(delassus_expression, enforce_symmetry);
    }

    /// \brief Constructor from a DelassusOperatorRigidBodySystemsTpl.
    /// \note The constructor does not compute the cholesky decomposition of the delassus.
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintModel,
      template<typename T> class StorageHolder>
    explicit DelassusOperatorDenseTpl(
      const DelassusOperatorRigidBodySystemsTpl<
        Scalar,
        Options,
        JointCollectionTpl,
        ConstraintModel,
        StorageHolder> & delassus_rigid_body,
      const bool enforce_symmetry = false)
    : DelassusOperatorDenseTpl()
    {
      rebuild(delassus_rigid_body, enforce_symmetry);
    }

    /// \brief Rebuild the internal data structure from a given matrix.
    /// \note This resets the physical compliance and numerical damping.
    template<typename MatrixDerived>
    void rebuild(const Eigen::MatrixBase<MatrixDerived> & mat)
    {
      PINOCCHIO_THROW_IF(mat.rows() != mat.cols(), std::runtime_error, "Expected a square matrix.");

      // copy input matrix
      m_delassus_matrix_storage.resize(mat.rows(), mat.cols());
      m_delassus_matrix = mat;

      // resize cholesky decomposition data
      m_cholesky_decomposition_data_storage.resize(mat.rows(), mat.cols());
      m_cholesky_decomposition_data.setZero();

      // resize/reset damping and compliance
      m_damping = BlockDiagonalMatrix::Zero(mat.rows());
      m_compliance_storage.resize(mat.rows());
      m_compliance.setZero();

      // mark decomposition as dirty
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Rebuild the internal data structure from a DelassusCholeskyExpressionTpl.
    template<typename ConstraintCholeskyDecomposition>
    void rebuild(
      const DelassusCholeskyExpressionTpl<ConstraintCholeskyDecomposition> & delassus_expression,
      const bool enforce_symmetry = false)
    {
      assert(delassus_expression.rows() == delassus_expression.size());
      assert(delassus_expression.cols() == delassus_expression.size());

      // retrieve delassus matrix
      const auto size = delassus_expression.size();
      m_delassus_matrix_storage.resize(size, size);
      delassus_expression.undampedMatrix(m_delassus_matrix, enforce_symmetry);
      m_delassus_matrix.diagonal() -= delassus_expression.getCompliance();

      // resize cholesky decomposition data
      m_cholesky_decomposition_data_storage.resize(size, size);
      m_cholesky_decomposition_data.setZero();

      // resize/reset damping and compliance
      m_damping = delassus_expression.getDamping();
      m_compliance_storage.resize(size);
      m_compliance = delassus_expression.getCompliance();

      // mark decomposition as dirty
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Rebuild the internal data structure from a DelassusOperatorRigidBodySystemsTpl.
    template<
      template<typename, int> class JointCollectionTpl,
      typename ConstraintModel,
      template<typename T> class StorageHolder>
    void rebuild(
      const DelassusOperatorRigidBodySystemsTpl<
        Scalar,
        Options,
        JointCollectionTpl,
        ConstraintModel,
        StorageHolder> & delassus_rigid_body,
      const bool enforce_symmetry = false)
    {
      assert(delassus_rigid_body.rows() == delassus_rigid_body.size());
      assert(delassus_rigid_body.cols() == delassus_rigid_body.size());

      // retrieve delassus matrix
      const auto size = delassus_rigid_body.size();
      m_delassus_matrix_storage.resize(size, size);
      delassus_rigid_body.undampedMatrix(m_delassus_matrix, enforce_symmetry);
      m_delassus_matrix.diagonal() -= delassus_rigid_body.getCompliance();

      // resize cholesky decomposition data
      m_cholesky_decomposition_data_storage.resize(size, size);
      m_cholesky_decomposition_data.setZero();

      // resize/reset damping and compliance
      m_damping = delassus_rigid_body.getDamping();
      m_compliance_storage.resize(size);
      m_compliance = delassus_rigid_body.getCompliance();

      // mark decomposition as dirty
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Comparison operator.
    bool operator==(const Self & other) const
    {
      if (&other == this)
        return true;

      return m_delassus_matrix == other.m_delassus_matrix                              //
             && m_cholesky_decomposition_data == other.m_cholesky_decomposition_data   //
             && m_cholesky_decomposition_dirty == other.m_cholesky_decomposition_dirty //
             && m_damping == other.m_damping                                           //
             && m_compliance == other.m_compliance;
    }

    /// \brief Comparison operator.
    bool operator!=(const Self & other) const
    {
      return !(*this == other);
    }

    /// \brief Evaluates the product delassus * x and stores it in res.
    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x,
      const Eigen::MatrixBase<MatrixOut> & res_,
      bool with_damping = true) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), size());
      MatrixOut & res = res_.const_cast_derived();
      res.noalias() = m_delassus_matrix * x;
      res.noalias() += m_compliance.asDiagonal() * x;
      if (with_damping)
      {
        m_damping.template applyOnTheRight<pinocchio::internal::add_assign_op>(x, res);
      }
    }

    /// \brief Update physical compliance from a vector.
    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance_vector)
    {
      m_compliance = compliance_vector;
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Update physical compliance from a scalar.
    void updateCompliance(const Scalar & compliance)
    {
      updateCompliance(Vector::Constant(size(), compliance));
    }

    /// \brief Update numerical damping from a vector.
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & damping_vector)
    {
      m_damping = damping_vector.asDiagonal();
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Update numerical damping from a scalar.
    void updateDamping(const Scalar & mu)
    {
      m_damping = std::move(BlockDiagonalMatrix::ScalarIdentity(size(), mu));
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Update numerical damping by copying an input block diagonal matrix.
    template<int OtherOptions, std::size_t OtherAlignment>
    void updateDamping(
      const BlockDiagonalMatrixTpl<Scalar, OtherOptions, OtherAlignment> &
        block_diagonal_damping_matrix)
    {
      if (&block_diagonal_damping_matrix == &m_damping)
        return;

      m_damping = block_diagonal_damping_matrix;
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Update numerical damping by moving an input block diagonal matrix.
    template<int OtherOptions, std::size_t OtherAlignment>
    void updateDamping(
      BlockDiagonalMatrixTpl<Scalar, OtherOptions, OtherAlignment> && block_diagonal_damping_matrix)
    {
      if (&block_diagonal_damping_matrix == &m_damping)
        return;

      m_damping = std::move(block_diagonal_damping_matrix);
      m_cholesky_decomposition_dirty = true;
    }

    /// \brief Updates the cholesky decomposition of the delassus.
    void updateDecomposition()
    {
      if (m_cholesky_decomposition_dirty)
      {
        m_cholesky_decomposition_data = m_delassus_matrix;
        m_damping.addTo(m_cholesky_decomposition_data);
        m_cholesky_decomposition_data += m_compliance.asDiagonal();
        computeCholeskyDecomposition();
        m_cholesky_decomposition_dirty = false;
      }
    }

    /// \brief Returns true if updateDecomposition() needs to be called in order to call
    /// solveInPlace.
    bool isDirty() const
    {
      return m_cholesky_decomposition_dirty;
    }

    /// \brief Fills input matrix with the dense representation of the Delassus.
    template<typename MatrixType>
    void matrix(
      const Eigen::MatrixBase<MatrixType> & mat,
      bool enforce_symmetry = false,
      bool with_damping = true) const
    {
      MatrixType & mat_ = mat.const_cast_derived();
      mat_ = m_delassus_matrix;
      if (with_damping)
      {
        m_damping.addTo(mat_);
      }
      mat_ += m_compliance.asDiagonal();
      if (enforce_symmetry)
      {
        enforceSymmetry(mat_);
      }
    }

    /// \brief Returns the dense representation of the Delassus.
    Matrix matrix(bool enforce_symmetry = false, bool with_damping = true) const
    {
      Matrix res(rows(), cols());
      matrix(res, enforce_symmetry, with_damping);
      return res;
    }

    /// \brief Fills input matrix with the dense representation of the Delassus.
    /// The numerical damping is NOT taken into account here.
    template<typename MatrixType>
    void
    undampedMatrix(const Eigen::MatrixBase<MatrixType> & mat, bool enforce_symmetry = false) const
    {
      matrix(mat, enforce_symmetry, false /*no damping*/);
    }

    /// \brief Returns the dense representation of the Delassus.
    /// The numerical damping is NOT taken into account here.
    Matrix undampedMatrix(bool enforce_symmetry = false) const
    {
      Matrix res(rows(), cols());
      undampedMatrix(res, enforce_symmetry);
      return res;
    }

    /// \brief solveInPlace operation returning the results of the inverse of the Delassus operator
    /// times the input matrix mat.
    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      PINOCCHIO_THROW_IF(
        isDirty(), std::logic_error,
        "The DelassusOperator has dirty quantities. Please call updateDecomposition() first.");
      m_cholesky_decomposition.solveInPlace(mat.const_cast_derived());
    }

    /// \brief Same as \ref solveInPlace but returns the result.
    template<typename MatrixLike>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike)
      solve(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) res(mat);
      solveInPlace(res);
      return res;
    }

    /// \brief Same as \ref solveInPlace but stores the result in res.
    template<typename MatrixDerivedIn, typename MatrixDerivedOut>
    void solve(
      const Eigen::MatrixBase<MatrixDerivedIn> & x,
      const Eigen::MatrixBase<MatrixDerivedOut> & res) const
    {
      res.const_cast_derived() = x;
      solveInPlace(res.const_cast_derived());
    }

    /// \brief Returns the inverse of the damped delassus matrix.
    Matrix inverse() const
    {
      Matrix res = Matrix::Identity(size(), size());
      solveInPlace(res);
      return res;
    }

    /// \brief Returns the number of rows/cols of the Delassus.
    /// The delassus represents a size() x size() linear operator.
    Eigen::Index size() const
    {
      return m_delassus_matrix.rows();
    }

    /// \brief Returns the number of rows of the Delassus.
    Eigen::Index rows() const
    {
      return m_delassus_matrix.rows();
    }

    /// \brief Returns the number of cols of the Delassus.
    Eigen::Index cols() const
    {
      return m_delassus_matrix.cols();
    }

    /// \brief Const getter for physical compliance.
    typename VectorStorage::ConstRefConstMapType & getCompliance() const
    {
      return m_compliance_storage.const_map();
    }

    /// \brief Const getter for numerical damping.
    const BlockDiagonalMatrix & getDamping() const
    {
      return m_damping;
    }

  protected:
    /// \brief Compute the cholesky decomposition of the matrix contained
    /// in m_cholesky_decomposition_data and stores it in place.
    void computeCholeskyDecomposition()
    {
      if (
        m_cholesky_decomposition.cols() != m_cholesky_decomposition_data.cols()
        || m_cholesky_decomposition.rows() != m_cholesky_decomposition_data.rows())
      {
        // if the decomposition does not point to the data with the right size,
        // we recreate it.
        m_cholesky_decomposition.~CholeskyDecomposition();
        new (&m_cholesky_decomposition) CholeskyDecomposition(m_cholesky_decomposition_data);
      }
      else
      {
        // otherwise we run the decomposition algorithm on the internal decomposition data.
        m_cholesky_decomposition.compute(m_cholesky_decomposition_data);
      }
    }

    /// \brief Storage for the delassus matrix.
    MatrixStorage m_delassus_matrix_storage;
    MatrixStorageRefMapType m_delassus_matrix = m_delassus_matrix_storage.map();

    /// \brief Data where the cholesky decomposition is stored.
    MatrixStorage m_cholesky_decomposition_data_storage;
    MatrixStorageRefMapType m_cholesky_decomposition_data =
      m_cholesky_decomposition_data_storage.map();

    /// \brief Cholesky decomposition algorithm.
    CholeskyDecomposition m_cholesky_decomposition;

    /// \brief Boolean to signal wether or not updateDecomposition() should be called.
    bool m_cholesky_decomposition_dirty;

    /// \brief Block diagonal numerical damping.
    BlockDiagonalMatrix m_damping;

    /// \brief Physical compliance.
    VectorStorage m_compliance_storage;
    VectorStorageRefMapType m_compliance = m_compliance_storage.map();
  }; // struct DelassusOperatorDenseTpl

} // namespace pinocchio
