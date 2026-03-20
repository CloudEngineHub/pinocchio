//
// Copyright (c) 2024-2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/algorithm/delassus-operator.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/algorithm/delassus-operator.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  namespace internal
  {

    template<typename Derived>
    struct SimplicialCholeskyWrapper : public Derived
    {
      typedef Eigen::SimplicialCholeskyBase<Derived> Base;

      using Base::derived;
      using Base::m_diag;
      using Base::m_info;
      using Base::m_matrix;
      using Base::m_P;
      using Base::m_Pinv;

      template<typename Rhs, typename Dest>
      void _solve_impl(const Eigen::MatrixBase<Rhs> & b, Eigen::MatrixBase<Dest> & dest) const
      {
        //    eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for
        //    solving, you must first call either compute() or symbolic()/numeric()");
        //    eigen_assert(m_matrix.rows()==b.rows());

        if (m_info != Eigen::Success)
          return;

        typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Rhs) PlainMatrix;
        typedef typename PlainMatrix::Scalar Scalar;

        typedef Eigen::Map<PlainMatrix, EIGEN_DEFAULT_ALIGN_BYTES> MapPlainMatrix;
        MapPlainMatrix tmp = MapPlainMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, b.rows(), b.cols()));
        if (m_P.size() > 0)
          tmp.noalias() = m_P * b;
        else
          tmp = b;

        if (m_matrix.nonZeros() > 0) // otherwise L==I
          derived().matrixL().solveInPlace(tmp);

        if (m_diag.size() > 0)
          tmp = m_diag.asDiagonal().inverse() * tmp;

        if (m_matrix.nonZeros() > 0) // otherwise U==I
          derived().matrixU().solveInPlace(tmp);

        if (m_P.size() > 0)
          dest.noalias() = m_Pinv * tmp;
      }

    }; // SimplicialCholeskyWrapper

    template<typename SparseCholeskySolver>
    struct getSparseCholeskySolverBase;

    template<typename SparseCholeskySolver> //, typename Base = typename SparseCholeskySolver::Base>
    struct SparseSolveInPlaceMethod;

#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
    template<typename MatrixType, int UpLo, SparseFactorization_t Solver, bool EnforceSquare>
    struct SparseSolveInPlaceMethod<Eigen::AccelerateImpl<MatrixType, UpLo, Solver, EnforceSquare>>
    {
      typedef Eigen::AccelerateImpl<MatrixType, UpLo, Solver, EnforceSquare> SparseCholeskySolver;

      template<typename Rhs, typename Dest>
      static void run(
        const SparseCholeskySolver & solver,
        const Eigen::MatrixBase<Rhs> & mat,
        const Eigen::MatrixBase<Dest> & dest)
      {
        dest.const_cast_derived() = solver.solve(mat.derived());
      }
    };
#endif

    template<typename SparseCholeskySolver>
    struct SparseSolveInPlaceMethod
    {
      template<typename Rhs, typename Dest>
      static void run(
        const SparseCholeskySolver & solver,
        const Eigen::MatrixBase<Rhs> & mat,
        const Eigen::MatrixBase<Dest> & dest)
      {
        static_assert(
          std::is_base_of_v<
            Eigen::SimplicialCholeskyBase<SparseCholeskySolver>, SparseCholeskySolver>,
          "The solver is not a base of SimplicialCholeskyBase.");
        typedef SimplicialCholeskyWrapper<SparseCholeskySolver> CholeskyWrapper;

        const CholeskyWrapper & wrapper = reinterpret_cast<const CholeskyWrapper &>(solver);
        wrapper._solve_impl(mat, dest.const_cast_derived());
      }
    };

  } // namespace internal

  template<typename _Scalar, int _Options, class _SparseCholeskyDecomposition>
  struct traits<DelassusOperatorSparseTpl<_Scalar, _Options, _SparseCholeskyDecomposition>>
  {
    typedef _SparseCholeskyDecomposition CholeskyDecomposition;
    typedef typename CholeskyDecomposition::MatrixType SparseMatrix;
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    static constexpr int RowsAtCompileTime = Eigen::Dynamic;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;

    typedef const Vector & getDampingReturnType;
  };

  template<typename _Scalar, int _Options, class SparseCholeskyDecomposition>
  struct DelassusOperatorSparseTpl
  : DelassusOperatorBase<DelassusOperatorSparseTpl<_Scalar, _Options, SparseCholeskyDecomposition>>
  {

    typedef DelassusOperatorSparseTpl Self;
    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;
    static constexpr int RowsAtCompileTime = traits<Self>::RowsAtCompileTime;

    typedef typename traits<Self>::SparseMatrix SparseMatrix;
    typedef typename traits<Self>::Vector Vector;
    typedef typename traits<Self>::Matrix Matrix;
    typedef SparseCholeskyDecomposition CholeskyDecomposition;
    typedef DelassusOperatorBase<Self> Base;

    template<typename MatrixDerived>
    explicit DelassusOperatorSparseTpl(const Eigen::SparseMatrixBase<MatrixDerived> & mat)
    : Base()
    , m_delassus_matrix(mat)
    , m_damped_delassus_matrix(mat)
    , m_cholsky_decomposition(mat)
    , m_cholesky_decomposition_dirty(true)
    , m_damping(Vector::Zero(mat.rows()))
    , m_compliance(Vector::Zero(mat.rows()))
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), mat.cols());
    }

    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance_vector)
    {
      for (Eigen::Index k = 0; k < size(); ++k)
      {
        m_damped_delassus_matrix.coeffRef(k, k) += -m_compliance[k] + compliance_vector[k];
      }
      m_compliance = compliance_vector;
      m_cholesky_decomposition_dirty = true;
    }

    void updateCompliance(const Scalar & m_compliance)
    {
      updateCompliance(Vector::Constant(size(), m_compliance));
    }

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & damping_vector)
    {
      for (Eigen::Index k = 0; k < size(); ++k)
      {
        m_damped_delassus_matrix.coeffRef(k, k) += -m_damping[k] + damping_vector[k];
      }
      m_damping = damping_vector;
      m_cholesky_decomposition_dirty = true;
    }

    void updateDamping(const Scalar & mu)
    {
      updateDamping(Vector::Constant(size(), mu));
    }

    void updateDecomposition()
    {
      PINOCCHIO_EIGEN_MALLOC_SAVE_STATUS();
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
      m_cholsky_decomposition.factorize(m_damped_delassus_matrix);
      PINOCCHIO_EIGEN_MALLOC_RESTORE_STATUS();
      m_cholesky_decomposition_dirty = false;
    }

    bool isDirty() const
    {
      return m_cholesky_decomposition_dirty;
    }

    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      PINOCCHIO_THROW_IF(
        isDirty(), std::logic_error,
        "The DelassusOperator has dirty quantities. Please call updateDecomposition() first.");
      internal::SparseSolveInPlaceMethod<CholeskyDecomposition>::run(
        m_cholsky_decomposition, mat.derived(), mat.const_cast_derived());
    }

    template<typename MatrixLike>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike)
      solve(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) res(mat);
      solveInPlace(res);
      return res;
    }

    template<typename MatrixDerivedIn, typename MatrixDerivedOut>
    void solve(
      const Eigen::MatrixBase<MatrixDerivedIn> & x,
      const Eigen::MatrixBase<MatrixDerivedOut> & res) const
    {
      res.const_cast_derived() = x;
      m_cholsky_decomposition._solve_impl(x, res.const_cast_derived());
    }

    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x,
      const Eigen::MatrixBase<MatrixOut> & res_,
      bool with_damping = true) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), size());
      MatrixOut & res = res_.const_cast_derived();
      res.noalias() = m_delassus_matrix * x;
      res.array() += m_compliance.array() * x.array();
      if (with_damping)
      {
        res.array() += m_damping.array() * x.array();
      }
    }

    Eigen::Index size() const
    {
      return m_delassus_matrix.rows();
    }
    Eigen::Index rows() const
    {
      return m_delassus_matrix.rows();
    }
    Eigen::Index cols() const
    {
      return m_delassus_matrix.cols();
    }

    SparseMatrix matrix(bool enforce_symmetry = false) const
    {
      m_damped_delassus_matrix = m_delassus_matrix;
      m_damped_delassus_matrix += m_compliance.asDiagonal();
      m_damped_delassus_matrix += m_damping.asDiagonal();
      if (enforce_symmetry)
      {
        // TODO: enforce symmetry for sparse matrices
        PINOCCHIO_THROW(
          std::invalid_argument, "enforceSymmetry not implemented for sparse matrices");
      }
      return m_damped_delassus_matrix;
    }

    const Vector & getCompliance() const
    {
      return m_compliance;
    }

    const Vector & getDamping() const
    {
      return m_damping;
    }

    SparseMatrix inverse() const
    {
      SparseMatrix identity_matrix(size(), size());
      identity_matrix.setIdentity();
      SparseMatrix res = m_cholsky_decomposition.solve(identity_matrix);
      return res;
    }

  protected:
    SparseMatrix m_delassus_matrix;
    mutable SparseMatrix m_damped_delassus_matrix;
    CholeskyDecomposition m_cholsky_decomposition;
    bool m_cholesky_decomposition_dirty;
    Vector m_damping;
    Vector m_compliance;

  }; // struct DelassusOperatorSparseTpl

} // namespace pinocchio
