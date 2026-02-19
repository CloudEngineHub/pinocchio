//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_math_block_diagonal_matrix_hpp__
#define __pinocchio_math_block_diagonal_matrix_hpp__

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/assign-operators.hpp"
#include "pinocchio/math/matrix-block-type.hpp"
#include "pinocchio/math/matrix-block-element.hpp"

#include "pinocchio/math/block-diagonal-matrix-base.hpp"
#include "pinocchio/math/block-diagonal-matrix-expression.hpp"
#include "pinocchio/math/block-diagonal-matrix-inverse.hpp"
#include "pinocchio/math/block-diagonal-matrix-sum.hpp"

#include "pinocchio/container/matrix-stack.hpp"
#include "pinocchio/utils/size-in-bytes.hpp"

#include "pinocchio/macros.hpp"

#include <limits>
#include <vector>

namespace pinocchio
{

  /// @brief A block-diagonal matrix with scalar type and options from the current context.
  typedef BlockDiagonalMatrixTpl<context::Scalar, context::Options> BlockDiagonalMatrix;

  /**
   * @ingroup pinocchio_math
   * @brief A memory-efficient representation of a block-diagonal matrix.
   *
   * @tparam _Scalar The scalar type of the matrix elements (e.g., double).
   * @tparam _Options The Eigen options for matrix storage (e.g., Eigen::RowMajor).
   * @tparam _Alignment The memory alignment for the data blocks.
   *
   * @details This class represents a block-diagonal matrix by storing only the non-zero blocks.
   *          It owns the memory for these blocks in a contiguous `MatrixStack` and uses a vector
   *          of `MatrixBlockElement` objects to describe the structure and provide views into the
   * data. This approach is highly efficient for storage and for performing matrix operations.
   */
  template<typename _Scalar, int _Options, std::size_t _Alignment>
  struct BlockDiagonalMatrixTpl
  : BlockDiagonalMatrixBase<BlockDiagonalMatrixTpl<_Scalar, _Options, _Alignment>>
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;
    static constexpr std::size_t Alignment = _Alignment;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
    typedef Eigen::Map<Matrix, to_eigen_alignment(Alignment)> MatrixMap;
    typedef Eigen::Map<const Matrix, to_eigen_alignment(Alignment)> ConstMatrixMap;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Map<Vector, to_eigen_alignment(Alignment)> VectorMap;
    typedef Eigen::Map<const Vector, to_eigen_alignment(Alignment)> ConstVectorMap;

    typedef MatrixStackTpl<Matrix, Alignment> MatrixStack;
    typedef MatrixBlockElementTpl<MatrixMap> MatrixBlockElement;
    typedef MatrixBlockElementTpl<ConstMatrixMap> ConstMatrixBlockElement;
    /**
     * @brief Defaut constructor.
     */
    BlockDiagonalMatrixTpl() {};

    /**
     * @brief Constructs a block-diagonal matrix from a given block pattern.
     * @param[in] block_pattern A vector of MatrixBlockElement describing each diagonal block in
     * order.
     */
    explicit BlockDiagonalMatrixTpl(const std::vector<MatrixBlockElement> & block_pattern);

    /**
     * @brief Constructs a block-diagonal matrix from an Eigen diagonal matrix expression.
     *
     * @tparam DiagonalVectorType The type of the underlying vector in the `Eigen::DiagonalWrapper`.
     *
     * @param[in] diagonal_expression An Eigen diagonal matrix expression, typically the result of
     *                                calling `.asDiagonal()` on an Eigen vector.
     *
     * @details This constructor converts a standard Eigen diagonal matrix into the specialized
     *          `BlockDiagonalMatrixTpl` representation. The resulting matrix will be composed of
     *          `N` blocks of size 1x1, where `N` is the dimension of the input matrix.
     *
     *          This is a deep-copy operation: the diagonal coefficients from the input
     *          expression are copied into the internal storage of this object.
     *
     * @code
     *   Eigen::Vector3d vec(1.0, 2.0, 3.0);
     *   // Create a BlockDiagonalMatrix from a standard Eigen diagonal matrix expression
     *   pinocchio::BlockDiagonalMatrix block_diag_matrix(vec.asDiagonal());
     *
     *   // block_diag_matrix now represents diag(1,2,3) using three 1x1 blocks.
     * @endcode
     */
    template<typename DiagonalVectorType>
    explicit BlockDiagonalMatrixTpl(
      const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression);

    /**
     * @brief Copy constructor.
     *
     * @param[in] other The BlockDiagonalMatrixTpl object to copy from.
     *
     * @details Performs a deep copy of the block-diagonal matrix. The new object will
     *          have its own independent copy of the underlying data in the `MatrixStack`.
     *          The `MatrixBlockElement` vector is also copied, and the internal `Eigen::Map`s
     *          are correctly re-pointed to view the newly allocated memory.
     *
     * @note This operation can be expensive as it involves allocating memory and copying
     *       all the data from the non-trivial matrix blocks.
     */
    BlockDiagonalMatrixTpl(const BlockDiagonalMatrixTpl & other)
    {
      *this = other;
    }

    /**
     * @brief Copy constructor from a block-diagonal matrix expression.
     *
     * @tparam Derived The derived type of the block-diagonal matrix expression.
     *
     * @param[in] other A block-diagonal matrix expression to assign from.
     *
     * @details This operator enables assignment from any type that inherits from
     *          `BlockDiagonalMatrixExpression`, such as `Inverse<BlockDiagonalMatrixTpl>`.
     *          The assignment is performed by calling `evalTo()` on the expression,
     *          which materializes the result into this matrix.
     *
     *          This enables lazy evaluation patterns where intermediate results
     *          (like matrix inverses) are not computed until they are assigned
     *          to a concrete storage type.
     *
     * @code
     *   pinocchio::BlockDiagonalMatrix A = ...;
     *   pinocchio::BlockDiagonalMatrix A_inv = A.inverse();
     * @endcode
     */
    template<typename Derived>
    BlockDiagonalMatrixTpl(const BlockDiagonalMatrixExpression<Derived> & other)
    {
      *this = other;
    }

    /**
     * @brief Move constructor.
     *
     * @param[in,out] other The BlockDiagonalMatrixTpl object to move from. After the move,
     *                      `other` is left in an invalid state.
     *
     * @details Transfers ownership of the underlying matrix data (`MatrixStack`) and block
     *          information from `other` to this object. This is a very efficient,
     *          constant-time operation as it avoids any memory allocation or data copying.
     */
    BlockDiagonalMatrixTpl(BlockDiagonalMatrixTpl && other) = default;

    /**
     * @brief Copy-assignment operator.
     *
     * @param[in] other The BlockDiagonalMatrixTpl object to copy from.
     * @return A reference to `*this` after the assignment.
     *
     * @details Replaces the contents of this instance with a deep copy of `other`.
     *          The existing data in `*this` is discarded. The new instance will have its
     *          own independent copy of the underlying data (`MatrixStack`) and block
     *          information. The internal `Eigen::Map`s are correctly re-pointed to
     *          view the newly allocated memory.
     *
     * @note This operation can be expensive if a memory reallocation is required.
     */
    BlockDiagonalMatrixTpl & operator=(const BlockDiagonalMatrixTpl & other);

    /**
     * @brief Move-assignment operator.
     *
     * @param[in,out] other The BlockDiagonalMatrixTpl object to move from. After the move,
     *                      `other` is left in a valid but unspecified state.
     * @return A reference to `*this` after the assignment.
     *
     * @details Transfers ownership of the underlying matrix data and block information
     *          from `other` to `*this`. The existing data in `*this` is properly released.
     *          This is a very efficient, constant-time operation that avoids any
     *          memory allocation or data copying.
     */
    BlockDiagonalMatrixTpl & operator=(BlockDiagonalMatrixTpl && other) = default;

    /**
     * @brief Assigns an Eigen diagonal matrix expression to this block-diagonal matrix.
     *
     * @tparam DiagonalVectorType The type of the underlying vector in the `Eigen::DiagonalWrapper`.
     *
     * @param[in] diagonal_expression An Eigen diagonal matrix expression, typically from
     * `.asDiagonal()`.
     *
     * @return A reference to `*this` after the assignment.
     *
     * @details This operator performs an **in-place value assignment**. It updates the numerical
     *          coefficients of the existing matrix blocks with the values from
     * `diagonal_expression`. It does **not** change the block structure (i.e., the number or sizes
     * of blocks) of this matrix.
     *
     *          **Preconditions:**
     *          - The dimensions of `diagonal_expression` must match the dimensions of `*this`.
     *          - The block structure of `*this` must be compatible with a diagonal matrix
     *            (i.e., it should be composed of 1x1 blocks).
     *
     * @warning This operator will result in undefined behavior if the dimensions do not match or if
     * the block structure is not diagonal.
     *
     * @code
     *   // Assume `block_diag_matrix` is already initialized, e.g., as a 3x3 matrix
     *   // of three 1x1 blocks.
     *   Eigen::Vector3d new_values(4, 5, 6);
     *   block_diag_matrix = new_values.asDiagonal(); // Updates the values in place
     * @endcode
     */
    template<typename DiagonalVectorType>
    BlockDiagonalMatrixTpl &
    operator=(const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression);

    /**
     * @brief Adds a diagonal matrix to this block-diagonal matrix in-place.
     * @tparam DiagonalVectorType The type of the underlying vector in the `Eigen::DiagonalWrapper`.
     * @param[in] diagonal_expression An Eigen diagonal matrix expression.
     * @return A reference to `*this` after the addition.
     */
    template<typename DiagonalVectorType>
    BlockDiagonalMatrixTpl &
    operator+=(const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression);

    /**
     * @brief Adds a diagonal matrix to this block-diagonal matrix.
     * @tparam DiagonalVectorType The type of the underlying vector in the `Eigen::DiagonalWrapper`.
     * @param[in] diagonal_expression An Eigen diagonal matrix expression.
     * @return A new BlockDiagonalMatrixTpl containing the sum.
     */
    template<typename DiagonalVectorType>
    Sum<BlockDiagonalMatrixTpl, Eigen::DiagonalWrapper<DiagonalVectorType>>
    operator+(const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression) const;

    /**
     * @brief Assignment operator from a block-diagonal matrix expression.
     *
     * @tparam Derived The derived type of the block-diagonal matrix expression.
     *
     * @param[in] other A block-diagonal matrix expression to assign from.
     *
     * @return A reference to `*this` after the assignment.
     *
     * @details This operator enables assignment from any type that inherits from
     *          `BlockDiagonalMatrixExpression`, such as `Inverse<BlockDiagonalMatrixTpl>`.
     *          The assignment is performed by calling `evalTo()` on the expression,
     *          which materializes the result into this matrix.
     *
     *          This enables lazy evaluation patterns where intermediate results
     *          (like matrix inverses) are not computed until they are assigned
     *          to a concrete storage type.
     *
     * @code
     *   pinocchio::BlockDiagonalMatrix A = ...;
     *   pinocchio::BlockDiagonalMatrix A_inv;
     *   A_inv = A.inverse();  // Computes and stores the inverse
     * @endcode
     */
    template<typename Derived>
    BlockDiagonalMatrixTpl & operator=(const BlockDiagonalMatrixExpression<Derived> & other)
    {
      other.evalTo(*this);
      return *this;
    }

    /**
     * @brief Checks for strict equality between two block-diagonal matrices.
     *
     * @param[in] other The other matrix to compare against.
     *
     * @return `true` if the matrices are equal, `false` otherwise.
     *
     * @details Two block-diagonal matrices are considered equal if and only if:
     *          1. Their overall dimensions (`rows` and `cols`) are identical.
     *          2. Their block patterns are identical (i.e., the `MatrixBlockElement` vectors are
     * the same).
     *          3. The numerical data in their underlying storage (`MatrixStack`) is
     * coefficient-wise equal.
     *
     * @note This comparison can be expensive as it may involve a full-data comparison of all
     * non-trivial blocks.
     */
    bool operator==(const BlockDiagonalMatrixTpl & other) const
    {
      // This implementation assumes MatrixBlockElementTpl has a valid operator==.
      // A correct implementation would need to handle the non-comparable Eigen::Map member.
      return m_rows == other.m_rows && m_cols == other.m_cols
             && m_matrix_block_elements == other.m_matrix_block_elements;
    }

    /**
     * @brief Checks for inequality between two block-diagonal matrices.
     *
     * @param[in] other The other matrix to compare against.
     *
     * @return `true` if the matrices are not equal, `false` otherwise.
     *
     * @details This operator is implemented as the negation of `operator==`.
     * @see operator==()
     */
    bool operator!=(const BlockDiagonalMatrixTpl & other) const
    {
      return !(*this == other);
    }

    /// @brief Checks if the matrix structure and its blocks are valid.
    bool isValid() const;

    /// @brief Returns the total number of rows of the full matrix.
    Eigen::Index rows() const
    {
      return m_rows;
    }

    /// @brief Returns the total number of columns of the full matrix.
    Eigen::Index cols() const
    {
      return m_cols;
    }

    /// @brief Returns the total number of elements in the full matrix (rows * cols).
    Eigen::Index size() const
    {
      return m_rows * m_cols;
    }

    /// @brief Returns the minimum coefficient of the matrix.
    Scalar minCoeff() const
    {
      PINOCCHIO_THROW_PRETTY_IF(
        blocks().empty(), std::runtime_error,
        "Unvalid use of minCoeff. You are using an empty block diagonal matrix");

      Scalar min_coeff = blocks()[0].map.minCoeff();
      for (const auto & block : blocks())
      {
        const Scalar block_min_coeff = block.map.minCoeff();
        if (block_min_coeff < min_coeff)
        {
          min_coeff = block_min_coeff;
        }
      }

      return min_coeff;
    }

    /// @brief Returns the maximum coefficient of the matrix.
    Scalar maxCoeff() const
    {
      PINOCCHIO_THROW_PRETTY_IF(
        blocks().empty(), std::runtime_error,
        "Unvalid use of maxCoeff. You are using an empty block diagonal matrix");

      Scalar max_coeff = blocks()[0].map.maxCoeff();
      for (const auto & block : blocks())
      {
        const Scalar block_max_coeff = block.map.maxCoeff();
        if (block_max_coeff > max_coeff)
        {
          max_coeff = block_max_coeff;
        }
      }

      return max_coeff;
    }

    /**
     * @brief Performs the matrix-matrix product `res = (*this) * rhs`.
     *
     * @tparam AssignOp  Assignment operation used to store the result
     *                   (e.g. ::pinocchio::internal::assign_op for direct assignment,
     *                   or ::pinocchio::internal::add_assign_op for accumulation).
     * @tparam MatrixRhs The Eigen type of the right-hand-side matrix.
     * @tparam MatrixRes The Eigen type of the result matrix.
     *
     * @param[in]  rhs The matrix to multiply with on the right.
     * @param[out] res The matrix where the result is stored. It must be pre-allocated
     *                 with dimensions `this->rows()` by `rhs.cols()`.
     *
     * @details This method computes the product of this block-diagonal matrix with a dense
     *          matrix `rhs`. It leverages the sparse structure of `*this` to perform the
     *          computation efficiently, avoiding unnecessary multiplications by zero.
     *
     *          As this version writes to a pre-allocated matrix, it avoids any dynamic
     *          memory allocation, making it suitable for real-time and performance-critical code.
     */
    template<
      typename AssignOp = pinocchio::internal::assign_op,
      typename MatrixRhs,
      typename MatrixRes>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixRhs> & rhs, const Eigen::MatrixBase<MatrixRes> & res) const;

    /**
     * @brief Performs the matrix-matrix product `res = lhs * (*this)`.
     *
     * @tparam AssignOp  Assignment operation used to store the result
     *                   (e.g. ::pinocchio::internal::assign_op for direct assignment,
     *                   or ::pinocchio::internal::add_assign_op for accumulation).
     * @tparam MatrixLhs The Eigen type of the left-hand-side matrix.
     * @tparam MatrixRes The Eigen type of the result matrix.
     *
     * @param[in]  lhs The matrix to multiply with on the left.
     * @param[out] res The matrix where the result is stored. It must be pre-allocated
     *                 with dimensions `this->cols()` by `lhs.rows()`.
     *
     * @details This method computes the product of this block-diagonal matrix with a dense
     *          matrix `lhs`. It leverages the sparse structure of `*this` to perform the
     *          computation efficiently, avoiding unnecessary multiplications by zero.
     *
     *          As this version writes to a pre-allocated matrix, it avoids any dynamic
     *          memory allocation, making it suitable for real-time and performance-critical code.
     */
    template<
      typename AssignOp = pinocchio::internal::assign_op,
      typename MatrixLhs,
      typename MatrixRes>
    void applyOnTheLeft(
      const Eigen::MatrixBase<MatrixLhs> & lhs, const Eigen::MatrixBase<MatrixRes> & res) const;

    /**
     * @brief Performs the matrix-matrix product `(*this) * rhs` and returns the result.
     *
     * @tparam MatrixRhs The Eigen type of the right-hand-side matrix.
     *
     * @param[in] rhs The matrix to multiply with on the right.
     *
     * @return A new matrix containing the result of the multiplication. The returned matrix
     *         type is deduced to be a plain, non-expression matrix.
     *
     * @details This is a convenience overload that allocates a new matrix to store the result
     *          of the product.
     *
     * @note For performance-critical applications, prefer the overload that accepts a
     *       pre-allocated result matrix to avoid repeated memory allocations.
     * @see void applyOnTheRight(const Eigen::MatrixBase<MatrixRhs> &, const
     * Eigen::MatrixBase<MatrixRes> &) const
     */
    template<typename MatrixRhs>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixRhs)
      applyOnTheRight(const Eigen::MatrixBase<MatrixRhs> & rhs) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixRhs) ReturnType;
      ReturnType res(rows(), rhs.cols()); // Assuming this should be rows(), rhs.cols()
      applyOnTheRight(rhs.derived(), res);
      return res;
    }

    /**
     * @brief Performs the matrix-matrix product `lhs * (*this)` and returns the result.
     *
     * @tparam MatrixLhs The Eigen type of the left-hand-side matrix.
     *
     * @param[in] lhs The matrix to multiply with on the left.
     *
     * @return A new matrix containing the result of the multiplication. The returned matrix
     *         type is deduced to be a plain, non-expression matrix.
     *
     * @details This is a convenience overload that allocates a new matrix to store the result
     *          of the product.
     *
     * @note For performance-critical applications, prefer the overload that accepts a
     *       pre-allocated result matrix to avoid repeated memory allocations.
     * @see void applyOnTheLeft(const Eigen::MatrixBase<MatrixLhs> &, const
     * Eigen::MatrixBase<MatrixRes> &) const
     */
    template<typename MatrixLhs>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLhs)
      applyOnTheLeft(const Eigen::MatrixBase<MatrixLhs> & lhs) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLhs) ReturnType;
      ReturnType res(lhs.rows(), cols()); // Assuming this should be rows(), rhs.cols()
      applyOnTheLeft(lhs.derived(), res);
      return res;
    }

    /**
     * @brief Performs the matrix-matrix product `(*this) * rhs` and returns the result.
     *
     * @tparam MatrixRhs The Eigen type of the right-hand-side matrix.
     *
     * @param[in] rhs The matrix to multiply with on the right.
     *
     * @return A new matrix containing the result of the multiplication. The returned matrix
     *         type is deduced to be a plain, non-expression matrix.
     */
    template<typename MatrixRhs>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixRhs)
    operator*(const Eigen::MatrixBase<MatrixRhs> & rhs) const
    {
      return applyOnTheRight(rhs);
    }

    /**
     * @brief Sets a dense matrix to be equal to this block-diagonal matrix.
     * @tparam Matrix An Eigen dense matrix type.
     * @param[out] matrix The dense matrix to be set. Its non-zero blocks will be filled, and
     *                    its off-diagonal blocks will be set to zero.
     */
    template<typename Matrix>
    void evalTo(const Eigen::MatrixBase<Matrix> & matrix) const;

    /**
     * @brief Adds this block-diagonal matrix to a dense matrix.
     * @tparam Matrix An Eigen dense matrix type.
     * @param[in,out] matrix The dense matrix to which this object will be added.
     */
    template<typename Matrix>
    void addTo(const Eigen::MatrixBase<Matrix> & matrix) const;

    /**
     * @brief Subtracts this block-diagonal matrix from a dense matrix.
     * @tparam Matrix An Eigen dense matrix type.
     * @param[in,out] matrix The dense matrix from which this object will be subtracted.
     */
    template<typename Matrix>
    void subTo(const Eigen::MatrixBase<Matrix> & matrix) const;

    /// @brief Returns a dense `Eigen::Matrix` representation of this block-diagonal matrix.
    /// @note This involves a memory allocation and data copy. For performance, prefer
    ///       operations like `evalTo` that work on existing memory.
    Matrix matrix() const;

    /**
     * @brief Fills a pre-allocated dense matrix with the values of this block-diagonal matrix.
     * @see evalTo
     * @tparam Matrix An Eigen dense matrix type.
     * @param[out] matrix The dense matrix to fill. Must have the correct dimensions.
     */
    template<typename Matrix>
    void matrix(const Eigen::MatrixBase<Matrix> & matrix) const;

    /// @brief Gets a const reference to the underlying memory stack containing block data.
    const MatrixStack & getMatrixStack() const
    {
      return m_matrix_stack;
    }

    /// @brief Gets a mutable reference to the underlying memory stack containing block data.
    MatrixStack & getMatrixStack()
    {
      return m_matrix_stack;
    }

    /// @brief Gets a const reference to the vector of block descriptors.
    const std::vector<MatrixBlockElement> & getMatrixBlockElements() const
    {
      return m_matrix_block_elements;
    }

    /// @brief Gets a mutable reference to the vector of block descriptors.
    std::vector<MatrixBlockElement> & getMatrixBlockElements()
    {
      return m_matrix_block_elements;
    }

    /// @copydoc getMatrixBlockElements
    std::vector<MatrixBlockElement> & blocks()
    {
      return m_matrix_block_elements;
    }

    /// @copydoc getMatrixBlockElements
    const std::vector<MatrixBlockElement> & blocks() const
    {
      return m_matrix_block_elements;
    }

    /**
     * @brief Fills a pre-allocated vector with the main diagonal of this block-diagonal matrix.
     *
     * @tparam DiagonalVector The Eigen type of the destination vector, slice, or expression.
     *
     * @param[out] diagonal_elements A pre-allocated, vector-like Eigen object that will be
     *                               filled with the diagonal elements. It **must** have a size
     *                               equal to the number of rows of this matrix.
     *
     * @details This is the core, high-performance method for extracting the full matrix diagonal.
     *          It avoids any memory allocation by writing directly into the provided destination.
     *
     *          The method iterates through the sequence of diagonal blocks. For each block,
     *          it computes its diagonal and copies it into the appropriate segment of the
     *          `diagonal_elements` vector.
     *
     * @see Vector diagonal() const
     */
    template<typename DiagonalVector>
    void diagonal(const Eigen::MatrixBase<DiagonalVector> & diagonal_elements) const;

    /**
     * @brief Extracts the main diagonal of this block-diagonal matrix into a new dense vector.
     *
     * @return A new dense column vector of size `rows()` containing the diagonal elements.
     *
     * @details This is a convenience method that allocates a new vector and calls the in-place
     *          `diagonal()` overload to fill it.
     *
     * @note For performance-critical code where repeated memory allocations should be avoided,
     *       prefer the overload that fills a pre-allocated vector.
     * @see void diagonal(const Eigen::MatrixBase<DiagonalVector>&) const
     */
    Vector diagonal() const
    {
      Vector diagonal_elements(rows());
      diagonal(diagonal_elements);
      return diagonal_elements;
    }

    /**
     * @brief Creates a square zero matrix represented as a `BlockDiagonalMatrixTpl`.
     *
     * @details This is a static factory method that constructs a `BlockDiagonalMatrixTpl` of
     *          dimensions `size` x `size` that represents a zero matrix.
     *
     *          The resulting object is represented with optimal efficiency, using a single
     *          block of type `BlockType::Zero` spanning the entire matrix dimension. This is
     *          a very lightweight operation as it does not require allocating memory for the
     *          matrix coefficients themselves.
     *
     * @param[in] size The dimension (number of rows and columns) for the resulting square
     *                 zero matrix. Must be non-negative.
     *
     * @return A `BlockDiagonalMatrixTpl` instance representing a `size` x `size` zero matrix.
     *
     * @code
     *   const Eigen::Index matrix_size = 10;
     *   auto zero_mat = pinocchio::BlockDiagonalMatrix::Zero(matrix_size);
     *
     *   // Check the properties of the created matrix
     *   // assert(zero_mat.rows() == matrix_size);
     *   // assert(zero_mat.cols() == matrix_size);
     *   // assert(zero_mat.matrix() == Eigen::MatrixXd::Zero(matrix_size, matrix_size));
     * @endcode
     */
    static BlockDiagonalMatrixTpl Zero(const Eigen::Index size);

    /**
     * @brief Creates a square matrix that is a scalar multiple of the identity matrix (s*I).
     *
     * @details This is a static factory method that constructs a `BlockDiagonalMatrixTpl` of
     *          dimensions `size` x `size` where all diagonal elements are equal to `value`
     *          and all off-diagonal elements are zero.
     *
     *          The resulting object is represented with optimal efficiency, using a single
     *          block of type `BlockType::Scalar`. This is highly memory-efficient as it only
     *          requires storing the single scalar `value`, regardless of the matrix size.
     *
     * @param[in] size The dimension (number of rows and columns) for the resulting square
     *                 matrix. Must be non-negative.
     * @param[in] value The scalar value to place on the main diagonal.
     *
     * @return A `BlockDiagonalMatrixTpl` instance representing a `size` x `size` scalar
     *         identity matrix.
     *
     * @code
     *   const Eigen::Index matrix_size = 10;
     *   const double scalar_value = 5.0;
     *   auto scalar_id_mat = pinocchio::BlockDiagonalMatrix::ScalarIdentity(matrix_size,
     * scalar_value);
     *
     *   // Check the properties of the created matrix
     *   // assert(scalar_id_mat.rows() == matrix_size);
     *   // assert(scalar_id_mat.cols() == matrix_size);
     *   // Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(matrix_size, matrix_size) *
     * scalar_value;
     *   // assert(scalar_id_mat.matrix().isApprox(expected));
     * @endcode
     */
    static BlockDiagonalMatrixTpl ScalarIdentity(const Eigen::Index size, const Scalar & value);

    /**
     * @brief Rebuilds the block-diagonal matrix from a given block pattern.
     * @param[in] new_block_pattern A vector of MatrixBlockElement describing each diagonal block in
     * order.
     */
    template<typename _MatrixBlockElement>
    void rebuild(const std::vector<_MatrixBlockElement> & new_block_pattern);

    /**
     * @brief Rebuilds the block-diagonal matrix from a given block pattern (pointer and size).
     * @param[in] new_block_pattern Pointer to an array of MatrixBlockElement.
     * @param[in] size Size of the array.
     */
    template<typename _MatrixBlockElement>
    void rebuild(const _MatrixBlockElement * new_block_pattern, const size_t size);

    /**
     * @brief Rebuilds the block-diagonal matrix from a diagonal expression.
     * @param[in] diagonal_expression An expression of the diagonal of the matrix.
     */
    template<typename DiagonalVectorType>
    void rebuild(const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression);

    /// \brief Returns a pointer to the underlying array serving as element storage.
    void * data()
    {
      return m_matrix_stack.data();
    }

    /// \brief Returns a pointer to the underlying array serving as element storage.
    const void * data() const
    {
      return m_matrix_stack.data();
    }

    /// \brief Returns the current memory footprint of this object in bytes.
    /// \details Sums up the sizes of all internal data members.
    std::size_t sizeInBytes() const
    {
      return 2 * ::pinocchio::sizeInBytes<Eigen::Index>()
             + m_matrix_stack
                 .sizeInBytes(); // TODO(jcarpent) complete + sizeInBytes(m_matrix_block_elements);
    }

    /// \brief Returns true if any coefficient (element) of this blocl element is NaN
    /// (Not‑a‑Number).
    bool hasNaN() const;

    /// \brief Returns an expression representing the inverse of this block diagonal matrix.
    Inverse<BlockDiagonalMatrixTpl> inverse() const;

  protected:
    /**
     * @brief Constructs a block-diagonal matrix from a given block pattern.
     * @param[in] block_pattern A vector of MatrixBlockElement describing each diagonal block in
     * order.
     */
    template<typename _MatrixBlockElement>
    void init_or_rebuild(const std::vector<_MatrixBlockElement> & block_pattern);

    template<typename _MatrixBlockElement>
    void init_or_rebuild(const _MatrixBlockElement * block_pattern, const size_t size);

    template<typename DiagonalVectorType>
    void init_or_rebuild(const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_expression);

    /**
     * @brief Clear the internal data structure for init or rebuild of the block diagonal matrix.
     */
    void clear();

    /**
     * @brief Generic implementation for assignment-like operations (e.g., =, +=, -=).
     * @tparam AssignOp Functor type that performs the assignment (e.g.,
     * `pinocchio::internal::SetTo`).
     * @tparam Matrix Dense Eigen matrix type.
     * @param[in,out] matrix The target matrix for the operation.
     */
    template<typename AssignOp, typename Matrix>
    void assign_op(const Eigen::MatrixBase<Matrix> & matrix) const;

  protected:
    /// @brief Total number of rows of the composite matrix.
    Eigen::Index m_rows = -1;

    /// @brief Total number of columns of the composite matrix.
    Eigen::Index m_cols = -1;

    /// @brief Contiguous memory storage for all non-trivial matrix blocks.
    MatrixStack m_matrix_stack;

    /// @brief A vector describing the sequence and properties of each diagonal block.
    std::vector<MatrixBlockElement> m_matrix_block_elements;
  }; // struct BlockDiagonalMatrixTpl

} // namespace pinocchio

#include "pinocchio/math/block-diagonal-matrix.hxx"

#endif // #ifndef __pinocchio_math_block_diagonal_matrix_hpp__
