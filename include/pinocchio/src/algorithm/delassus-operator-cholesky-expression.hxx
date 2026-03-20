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
  // Forward declaration of Unsafe specialization.
  template<typename _ContactCholeskyDecomposition>
  struct Unsafe<DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition>>;

  template<typename ContactCholeskyDecomposition>
  struct traits<DelassusCholeskyExpressionTpl<ContactCholeskyDecomposition>>
  {
    static constexpr int RowsAtCompileTime = Eigen::Dynamic;
    typedef typename ContactCholeskyDecomposition::Scalar Scalar;
    typedef typename ContactCholeskyDecomposition::Matrix Matrix;
    typedef typename ContactCholeskyDecomposition::Vector Vector;

    typedef typename ContactCholeskyDecomposition::EigenStorageVector EigenStorageVector;
    typedef typename ContactCholeskyDecomposition::BlockDiagonalMatrix BlockDiagonalMatrix;
    typedef const BlockDiagonalMatrix & getDampingReturnType;
  };

  /// \brief Unsafe version of DelassusCholeskyExpressionTpl.
  /// Allows direct access to protected members for expert users.
  template<typename _ContactCholeskyDecomposition>
  struct Unsafe<DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition>>
  {
    typedef DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition> SafeSelf;
    typedef typename SafeSelf::BlockDiagonalMatrix BlockDiagonalMatrix;

    explicit Unsafe(SafeSelf & self)
    : self(self)
    {
    }

    /// \brief Signal the delassus that updateDecomposition() should be called.
    /// This is typically called after damping has been directly modified via damping().
    void makeDirty()
    {
      self.self.updateSumComplianceDamping();
    }

    /// \brief Direct access to the block diagonal damping.
    BlockDiagonalMatrix & damping()
    {
      return self.self.m_damping;
    }

  protected:
    SafeSelf & self;
  };

  // TODO(jcarpent): change const_cast usage.
  template<typename _ContactCholeskyDecomposition>
  struct DelassusCholeskyExpressionTpl
  : DelassusOperatorBase<DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition>>
  {
    typedef _ContactCholeskyDecomposition ContactCholeskyDecomposition;
    typedef typename ContactCholeskyDecomposition::Scalar Scalar;
    typedef typename ContactCholeskyDecomposition::Vector Vector;
    typedef typename ContactCholeskyDecomposition::Matrix Matrix;
    typedef typename ContactCholeskyDecomposition::RowMatrix RowMatrix;
    typedef DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition> Self;
    typedef DelassusOperatorBase<Self> Base;
    typedef typename ContactCholeskyDecomposition::EigenStorageVector EigenStorageVector;
    typedef typename ContactCholeskyDecomposition::BlockDiagonalMatrix BlockDiagonalMatrix;
    static constexpr int Options = ContactCholeskyDecomposition::Options;
    typedef DelassusOperatorDenseTpl<Scalar, Options> DelassusOperatorDense;

    typedef
      typename SizeDepType<Eigen::Dynamic>::template BlockReturn<RowMatrix>::Type RowMatrixBlockXpr;
    typedef typename SizeDepType<Eigen::Dynamic>::template BlockReturn<RowMatrix>::ConstType
      RowMatrixConstBlockXpr;

    static constexpr int RowsAtCompileTime =
      traits<DelassusCholeskyExpressionTpl>::RowsAtCompileTime;

    /// \brief Cast this class to its unsafe version.
    Unsafe<Self> unsafe()
    {
      return Unsafe<Self>(*this);
    }
    friend struct Unsafe<Self>;

    /// \brief Default constructor from a cholesky decomposition.
    explicit DelassusCholeskyExpressionTpl(ContactCholeskyDecomposition & self)
    : Base()
    , self(self)
    {
    }

    /// \brief Evaluates the product Delassus * x and stores it in res.
    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x,
      const Eigen::MatrixBase<MatrixOut> & res,
      bool with_damping = true) const
    {

      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), self.constraintDim());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), self.constraintDim());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), x.cols());

      res.const_cast_derived().noalias() = self.delassus_block * x;
      if (with_damping)
      {
        self.m_sum_compliance_damping.template applyOnTheRight<pinocchio::internal::add_assign_op>(
          x, res.const_cast_derived());
      }
      else
      {
        // take only compliance into account
        res.const_cast_derived().noalias() += self.compliance.asDiagonal() * x;
      }

      // const auto U1 = self.U.topLeftCorner(self.constraintDim(), self.constraintDim());
      // {
      //   typedef Eigen::Map<RowMatrix> MapType;
      //   MapType tmp_mat = MapType(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, x.rows(), x.cols()));
      //   //            tmp_mat.noalias() = U1.adjoint() * x;
      //   triangularMatrixMatrixProduct<Eigen::UnitLower>(U1.adjoint(), x.derived(), tmp_mat);

      //   // The following commented lines produced some memory allocation.
      //   // Should be replaced by a manual loop
      //   //          tmp_mat.array().colwise() *= -self.D.head(self.constraintDim()).array();
      //   for (Eigen::Index i = 0; i < x.cols(); ++i)
      //     tmp_mat.col(i).array() *= -self.D.head(self.constraintDim()).array();

      //   //            res.const_cast_derived().noalias() = U1 * tmp_mat;
      //   triangularMatrixMatrixProduct<Eigen::UnitUpper>(U1, tmp_mat, res.const_cast_derived());
      // }
    }

    ///
    /// \brief Update the decomposition after a call to updateDamping or updateCompliance.
    ///
    /// \remarks isDirty() allows to retrieve the current status of the decomposition.
    ///
    void updateDecomposition()
    {
      self.computeDelassusCholeskyDecomposition();
    }

    /// \brief Returns true if updateDecomposition() needs to be called in order to call
    /// solveInPlace.
    bool isDirty() const
    {
      return self.isDirty();
    }

    /// \brief solveInPlace operation returning the resultsof the inverse of the Delassus operator
    /// times the input x.
    template<typename MatrixDerived>
    void solveInPlace(const Eigen::MatrixBase<MatrixDerived> & x) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), self.constraintDim());

      PINOCCHIO_THROW_IF(
        self.isDirty(), std::logic_error,
        "The DelassusOperator has dirty quantities. Please call updateDecomposition() first.");

      const auto U1 = self.U.topLeftCorner(self.constraintDim(), self.constraintDim())
                        .template triangularView<Eigen::UnitUpper>();

      U1.solveInPlace(x.const_cast_derived());

      // The following commented lines produced some memory allocation.
      // Should be replaced by a manual loop
      //        x.const_cast_derived().array().colwise() *=
      //        -self.Dinv.head(self.constraintDim()).array();
      for (Eigen::Index i = 0; i < x.cols(); ++i)
        x.const_cast_derived().col(i).array() *= -self.Dinv.head(self.constraintDim()).array();

      U1.adjoint().solveInPlace(x);
    }

    /// \brief Same as solveInPlace but stores the result in res.
    template<typename MatrixDerivedIn, typename MatrixDerivedOut>
    void solve(
      const Eigen::MatrixBase<MatrixDerivedIn> & x,
      const Eigen::MatrixBase<MatrixDerivedOut> & res) const
    {
      res.const_cast_derived() = x;
      solveInPlace(res.const_cast_derived());
    }

    template<typename MatrixDerived>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixDerived)
      solve(const Eigen::MatrixBase<MatrixDerived> & x) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixDerived) ReturnType;
      ReturnType res(self.constraintDim(), x.cols());
      solve(x.derived(), res);
      return res;
    }

    /// \brief Returns the Constraint Cholesky decomposition associated to this
    /// DelassusCholeskyExpression.
    const ContactCholeskyDecomposition & cholesky() const
    {
      return self;
    }

    /// \brief Returns the Constraint Cholesky decomposition associated to this
    /// DelassusCholeskyExpression.
    ContactCholeskyDecomposition & cholesky()
    {
      return self;
    }

    /// \brief Returns the matrix resulting from the decomposition.
    Matrix matrix(bool enforce_symmetry = false, bool with_damping = true) const
    {
      Matrix res = self.getInverseOperationalSpaceInertiaMatrix(enforce_symmetry);
      if (!with_damping)
      {
        getDamping().subTo(res);
      }
      return res;
    }

    /// \brief Fill the input matrix with the matrix resulting from the decomposition.
    template<typename MatrixType>
    void matrix(
      const Eigen::MatrixBase<MatrixType> & mat,
      bool enforce_symmetry = false,
      bool with_damping = true) const
    {
      self.getInverseOperationalSpaceInertiaMatrix(mat.const_cast_derived(), enforce_symmetry);
      if (!with_damping)
      {
        getDamping().subTo(mat.const_cast_derived());
      }
    }

    /// \brief Returns the matrix resulting from the decomposition.
    /// The numerical damping is NOT taken into account here.
    Matrix undampedMatrix(bool enforce_symmetry = false) const
    {
      Matrix res = matrix(enforce_symmetry, false /*no damping*/);
      return res;
    }

    /// \brief Fill the input matrix with the matrix resulting from the decomposition
    /// The numerical damping is NOT taken into account here.
    template<typename MatrixType>
    void
    undampedMatrix(const Eigen::MatrixBase<MatrixType> & mat, bool enforce_symmetry = false) const
    {
      matrix(mat, enforce_symmetry, false /*no damping*/);
    }

    ///
    /// \brief Returns the current compliance vector.
    ///
    const typename EigenStorageVector::ConstMapType getCompliance() const
    {
      return self.getCompliance();
    }

    ///
    /// \brief Returns the current damping as a block diagonal matrix.
    ///
    const BlockDiagonalMatrix & getDamping() const
    {
      return self.getDamping();
    }

    Matrix inverse() const
    {
      return self.getOperationalSpaceInertiaMatrix();
    }

    ///
    /// \brief Returns the corresponding dense delassus operator.
    ///
    DelassusOperatorDense dense(bool enforce_symmetry = false) const
    {
      return DelassusOperatorDense(*this, enforce_symmetry);
    }

    ///
    /// \brief Add a compliance term to the diagonal of the Delassus matrix. The compliance terms
    /// should be all positives.
    ///
    /// \param[in] compliances Vector of compliances related to the constraints.
    ///
    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliances)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateCompliance(compliances);
    }

    ///
    /// \brief Add a compliance term to the diagonal of the Delassus matrix. The compliance term
    /// should be positive.
    ///
    /// \param[in] compliance Compliance of the constraints.
    ///
    void updateCompliance(const Scalar & compliance)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateCompliance(compliance);
    }

    ///
    /// \brief Add a damping term to the diagonal of the Delassus matrix. The damping terms should
    /// be all positives.
    ///
    /// \param[in] mus Vector of positive regularization factor allowing to enforce the definite
    /// positiveness of the matrix.
    ///
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & mus)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateDamping(mus);
    }

    ///
    /// \brief Add a damping term to the diagonal of the Delassus matrix. The damping term should be
    /// positive.
    ///
    /// \param[in] mu Regularization factor allowing to enforce the definite positiveness of the
    /// matrix.
    ///
    void updateDamping(const Scalar & mu)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateDamping(mu);
    }

    ///
    /// \brief Update the damping from a block diagonal matrix (copy overload).
    ///
    template<int OtherOptions, std::size_t OtherAlignment>
    void updateDamping(
      const BlockDiagonalMatrixTpl<Scalar, OtherOptions, OtherAlignment> & block_damping)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateDamping(block_damping);
    }

    ///
    /// \brief Update the damping from a block diagonal matrix (move overload).
    ///
    template<int OtherOptions, std::size_t OtherAlignment>
    void
    updateDamping(BlockDiagonalMatrixTpl<Scalar, OtherOptions, OtherAlignment> && block_damping)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateDamping(std::move(block_damping));
    }

    /// \brief Returns the number of rows/cols of the Delassus.
    /// The delassus represents a size() x size() linear operator.
    Eigen::Index size() const
    {
      return self.constraintDim();
    }

    /// \brief Returns the number of rows of the Delassus.
    Eigen::Index rows() const
    {
      return size();
    }

    /// \brief Returns the number of cols of the Delassus.
    Eigen::Index cols() const
    {
      return size();
    }

    /// \brief Returns the current memory footprint of this object in bytes.
    /// \details Sums up the sizes of all internal data members.
    std::size_t sizeInBytes() const
    {
      return self.sizeInBytes();
    }

  protected:
    ContactCholeskyDecomposition & self;
  }; // DelassusCholeskyExpression

} // namespace pinocchio
