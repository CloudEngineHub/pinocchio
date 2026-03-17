//
// Copyright (c) 2026 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/math.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/math.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  // Forward declaration
  template<typename MatrixOrMap, typename Enable>
  struct MatrixBlockElementTpl;

  template<typename MapType>
  struct traits<MatrixBlockElementTpl<MapType>, std::enable_if_t<helper::is_eigen_map_v<MapType>>>
  {
    typedef typename MapType::PlainObject Matrix;
    typedef typename Matrix::Scalar Scalar;
    static constexpr auto Options = Matrix::Options;
    typedef MapType MatrixContainer;
    typedef Eigen::Matrix<Scalar, Matrix::RowsAtCompileTime, 1, Options> Vector;
    typedef MatrixBlockElementTpl<Matrix> PlainBlockElement;
  };

  template<typename MatrixType>
  struct traits<
    MatrixBlockElementTpl<MatrixType>,
    std::enable_if_t<helper::is_eigen_matrix_v<MatrixType>>>
  {
    typedef MatrixType Matrix;
    typedef typename Matrix::Scalar Scalar;
    static constexpr auto Options = Matrix::Options;
    typedef MatrixType MatrixContainer;
    typedef Eigen::Matrix<Scalar, Matrix::RowsAtCompileTime, 1, Options> Vector;
    typedef MatrixBlockElementTpl<Matrix> PlainBlockElement;
  };

  /**
   * @ingroup pinocchio_math
   * @brief A descriptor for a non-owning matrix block, implemented as a view using `Eigen::Map`.
   *
   * @tparam MapType The container type for the block. This specialization is enabled only when
   *                 `MapType` is an `Eigen::Map` type.
   *
   * @details This struct is a template specialization that describes a single block within a
   *          larger structured matrix. Its defining feature is that it does **not own** any
   *          memory for the matrix coefficients. Instead, it uses an `Eigen::Map` (`map` member)
   *          to provide a non-owning, mutable view into a region of memory that is managed
   *          externally (e.g., by a `MatrixStack` in a `BlockDiagonalMatrixTpl`).
   *
   *          It inherits metadata like `type` and `size` from the `MatrixBlockElementPlain` class.
   *
   * @note The copy constructor and copy-assignment operator perform a **shallow copy**. This means
   *       that both the original and the copied object will have their `map` members pointing
   *       to the **exact same underlying memory buffer**.
   */
  template<typename MapType>
  struct MatrixBlockElementTpl<MapType, std::enable_if_t<helper::is_eigen_map_v<MapType>>>
  : MatrixBlockElementPlain<MatrixBlockElementTpl<MapType>>
  {
    /// @brief The type of this specialized class.
    typedef MatrixBlockElementTpl Self;
    /// @brief The scalar type of the matrix elements (e.g., double).
    typedef typename traits<Self>::Scalar Scalar;
    /// @brief The non-owning container type, which is an `Eigen::Map`.
    typedef typename traits<Self>::MatrixContainer MatrixMap;
    /// @brief The equivalent dense Eigen matrix type.
    typedef typename traits<Self>::Matrix Matrix;
    /// @brief The equivalent dense Eigen vector type.
    typedef typename traits<Self>::Vector Vector;

    /// @brief The base class from which this class inherits.
    typedef MatrixBlockElementPlain<MatrixBlockElementTpl> Base;

    using Base::size;
    using Base::type;
    using Base::operator=;

    /// @brief An Eigen::Map that provides a non-owning view of the memory corresponding to this
    /// block.
    /// @note For block types that don't store data (e.g., `Identity`, `Zero`), this map will be
    /// null.
    MatrixMap map;

    /// @brief Default constructor. Initializes to an invalid state (Undefined type, size -1, null
    /// map).
    MatrixBlockElementTpl()
    : Base()
    , map(nullptr, 0, 0)
    {
    }

    /**
     * @brief Constructs a block descriptor for types that do not require external data.
     * @details This is used for structural blocks like `Identity` or `Zero`. The internal map
     *          remains null.
     * @param[in] type The structural type of the block.
     * @param[in] size The dimension of the (square) block.
     */
    MatrixBlockElementTpl(const MatrixBlockType type, const Eigen::Index size)
    : Base(type, size)
    , map(nullptr, 0, 0)
    {
    }

    /**
     * @brief Constructs a block descriptor for types that map to existing data.
     * @details This is used for data-backed blocks like `Plain`, `Diagonal`, or `ScalarIdentity`.
     * @param[in] type The structural type of the block.
     * @param[in] size The dimension of the (square) block.
     * @param[in] matrix_map An Eigen::Map providing a non-owning view of the block's data.
     */
    MatrixBlockElementTpl(
      const MatrixBlockType type, const Eigen::Index size, const MatrixMap matrix_map)
    : Base(type, size)
    , map(matrix_map)
    {
    }

    /**
     * @brief Default copy constructor (shallow copy).
     * @details Creates a copy of the block descriptor. The new object's `map` will view the
     *          **same memory** as the original. No matrix data is duplicated.
     */
    MatrixBlockElementTpl(const MatrixBlockElementTpl & other)
    : MatrixBlockElementTpl()
    {
      *this = other;
    }

    /**
     * @brief Default copy-assignment operator (shallow copy).
     * @details Assigns from another block descriptor. After assignment, this object's `map`
     *          will view the **same memory** as the other object.
     */
    MatrixBlockElementTpl & operator=(const MatrixBlockElementTpl & other)
    {
      if (this != &other)
      {
        Base::operator=(other);
        // For Eigen::Map, operator= only works if the copied map
        // has the same size than the current map.
        // If not, Eigen triggers an assert and fails at runtime.
        // Since the map we want to copy is not necessarily the same size
        // as the one we have, we recreate the map in place so that it
        // points to the same data as the other map.
        new (&map) MapType(other.map);
      }
      return *this;
    }

    /**
     * @brief Checks for strict equality between two block descriptors.
     *
     * @param[in] other The other block descriptor to compare against.
     * @return `true` if the blocks are strictly equal, `false` otherwise.
     *
     * @details Two block descriptors are considered equal if they have the same `type` and `size`,
     *          and the data viewed by their internal `map`s is **coefficient-wise equal**.
     *
     * @note This performs a deep, numerical comparison of the underlying data via
     * `Eigen::Map::operator==`, which can be computationally expensive. It does not simply compare
     * pointers.
     */
    bool operator==(const MatrixBlockElementTpl & other) const
    {
      if (this == &other)
        return true;
      return Base::operator==(other) && map == other.map;
    }

    /**
     * @brief Checks for inequality between two block descriptors.
     * @param[in] other The other block descriptor to compare against.
     * @return `true` if the blocks are not strictly equal, `false` otherwise.
     * @details This is the logical negation of `operator==`.
     * @see operator==()
     */
    bool operator!=(const MatrixBlockElementTpl & other) const
    {
      return !(*this == other);
    }

    /**
     * @brief Re-points the internal Eigen::Map to a new memory location.
     *
     * @tparam OtherMatrix The matrix type of the new map.
     * @tparam OtherAlignment The alignment option of the new map.
     * @tparam OtherStrideType The stride type of the new map.
     *
     * @param[in] other_map The new `Eigen::Map` to view.
     *
     * @details This is an advanced operation that uses placement-new to reconstruct the internal
     *          `map` member in-place, making it view the memory provided by `other_map`.
     */
    template<typename OtherMatrix, int OtherAlignment, typename OtherStrideType>
    void remap(Eigen::Map<OtherMatrix, OtherAlignment, OtherStrideType> & other_map)
    {
      new (&map) MatrixMap(other_map.data(), other_map.rows(), other_map.cols());
    }

    /**
     * @brief Checks if the block descriptor is valid and self-consistent.
     *
     * @details A block is considered valid if its base is valid (`size > 0` and `type !=
     * Undefined`) and its `map` points to valid data if the `type` requires it (e.g., `Plain`,
     * `Diagonal`).
     *
     * @return `true` if the block info is valid, `false` otherwise.
     */
    bool isValid() const
    {
      bool is_invalid = !Base::isValid() || (isDataBlock(type()) && map.data() == nullptr);
      return !is_invalid;
    }

    /// @brief Returns a mutable reference to the underlying Eigen::Map.
    MatrixMap & container()
    {
      return map;
    }

    /// @brief Returns a const reference to the underlying Eigen::Map.
    const MatrixMap & container() const
    {
      return map;
    }

    /**
     * @brief Gets a const pointer to the beginning of the map's data buffer.
     *
     * @return A const pointer to the first element (`Scalar*`) of the map's raw data.
     *
     * @details This provides direct, low-level read-only access to the matrix's coefficients.
     *          The data is stored in a contiguous block, typically in column-major order for
     *          Eigen matrices. This is useful for interoperability with C-style APIs or other
     *          libraries that operate on raw memory buffers.
     *
     * @warning The returned pointer is only valid as long as this object exists and its
     *          `m_map` member is not reallocated. Accessing it after the object is
     *          destroyed leads to undefined behavior.
     */
    const Scalar * data() const
    {
      return map.data();
    }

    /**
     * @brief Gets a mutable pointer to the beginning of the map's data buffer.
     * @copydoc data() const
     * @return A mutable pointer to the first element (`Scalar*`) of the map's raw data.
     */
    Scalar * data()
    {
      return map.data();
    }

  }; // struct MatrixBlockElementTpl

  /**
   * @ingroup pinocchio_math
   * @brief A descriptor for an owning matrix block, which stores its data in an Eigen::Matrix.
   *
   * @tparam MatrixType The container type for the block. This specialization is enabled only when
   *                    `MatrixType` is an owning `Eigen::Matrix` type (e.g., `Eigen::MatrixXd`).
   *
   * @details This struct is a template specialization that describes and **owns** a single block
   *          within a larger structured matrix. Unlike the `Eigen::Map` specialization which only
   *          provides a view, this version contains an `Eigen::Matrix` member (`m_matrix`) that
   *          holds the numerical coefficients of the block.
   *
   *          It inherits metadata like `type` and `size` from the `MatrixBlockElementPlain` class.
   *
   * @note The copy constructor and copy-assignment operator perform a **deep copy**. This means
   *       that the `m_matrix` member is fully duplicated, ensuring that the new object has its
   *       own independent copy of the data.
   */
  template<typename MatrixType>
  struct MatrixBlockElementTpl<MatrixType, std::enable_if_t<helper::is_eigen_matrix_v<MatrixType>>>
  : MatrixBlockElementPlain<MatrixBlockElementTpl<MatrixType>>
  {
    /// @brief The type of this specialized class.
    typedef MatrixBlockElementTpl Self;
    /// @brief The scalar type of the matrix elements (e.g., double).
    typedef typename traits<Self>::Scalar Scalar;
    /// @brief The owning container type, which is an `Eigen::Matrix`.
    typedef typename traits<Self>::MatrixContainer Matrix;
    // Note: The MatrixMap typedef from the other specialization is not applicable here.
    // We will assume `MatrixContainer` is `MatrixType`.
    // typedef typename traits<Self>::MatrixContainer MatrixMap;
    /// @brief The equivalent dense Eigen vector type.
    typedef typename traits<Self>::Vector Vector;

    /// @brief The base class from which this class inherits.
    typedef MatrixBlockElementPlain<MatrixBlockElementTpl> Base;

    using Base::size;
    using Base::type;
    using Base::operator=;

  protected:
    /// @brief The owning container for the matrix block's data.
    MatrixType m_matrix;

  public:
    /// @brief Default constructor. Initializes to an invalid state (Undefined type, size -1).
    MatrixBlockElementTpl()
    : Base()
    {
    }

    /**
     * @brief Constructs a block descriptor for types that may not require explicit data storage.
     * @details This can be used for structural blocks like `Identity` or `Zero`. The internal
     *          `m_matrix` member will be default-initialized.
     * @param[in] type The structural type of the block.
     * @param[in] size The dimension of the (square) block.
     */
    MatrixBlockElementTpl(const MatrixBlockType type, const Eigen::Index size)
    : Base(type, size)
    {
      switch (type)
      {
      case MatrixBlockType::Zero:
      case MatrixBlockType::Identity:
        break;
      case MatrixBlockType::ScalarIdentity:
        m_matrix.resize(1, 1);
        break;
      case MatrixBlockType::Diagonal:
        m_matrix.resize(size, 1);
        break;
      case MatrixBlockType::Plain:
        m_matrix.resize(size, size);
        break;
      default:
        assert(false && "Should never happened");
      }
    }

    /**
     * @brief Constructs a block descriptor and initializes its data by copying from another matrix.
     * @details This is the primary constructor for data-backed blocks like `Plain` or `Diagonal`.
     *          It performs a **deep copy** of the input matrix data into its internal storage.
     * @param[in] type The structural type of the block.
     * @param[in] size The dimension of the (square) block.
     * @param[in] matrix_data An Eigen matrix or expression whose data will be copied.
     */
    template<typename OtherMatrix>
    MatrixBlockElementTpl(
      const MatrixBlockType type,
      const Eigen::Index size,
      const Eigen::MatrixBase<OtherMatrix> & matrix_data)
    : Base(type, size)
    , m_matrix(matrix_data)
    {
    }

    /**
     * @brief Copy constructor (deep copy).
     * @details Creates a complete, independent copy of the block descriptor, including a full
     *          duplication of the underlying matrix data in `m_matrix`.
     */
    MatrixBlockElementTpl(const MatrixBlockElementTpl & other) = default;

    /**
     * @brief Copy-assignment operator (deep copy).
     * @details Replaces the content of this block with a full copy of the other block, including
     *          its matrix data.
     */
    MatrixBlockElementTpl & operator=(const MatrixBlockElementTpl & other) = default;

    /**
     * @brief Checks for strict equality between two block descriptors.
     *
     * @param[in] other The other block descriptor to compare against.
     * @return `true` if the blocks are strictly equal, `false` otherwise.
     *
     * @details Two block descriptors are considered equal if they have the same `type` and `size`,
     *          and their internal `m_matrix` members are **coefficient-wise equal**.
     */
    bool operator==(const MatrixBlockElementTpl & other) const
    {
      if (this == &other)
        return true;
      return Base::operator==(other) && m_matrix == other.m_matrix;
    }

    /**
     * @brief Checks for inequality between two block descriptors.
     * @param[in] other The other block descriptor to compare against.
     * @return `true` if the blocks are not strictly equal, `false` otherwise.
     * @details This is the logical negation of `operator==`.
     * @see operator==()
     */
    bool operator!=(const MatrixBlockElementTpl & other) const
    {
      return !(*this == other);
    }

    /**
     * @brief Checks if the block descriptor is valid.
     * @details A block is considered valid if its base is valid (`size > 0` and `type !=
     * Undefined`). Unlike the non-owning version, it doesn't need to check for null data pointers.
     * @return `true` if the block info is valid, `false` otherwise.
     */
    bool isValid() const
    {
      return Base::isValid();
    }

    /// @brief Returns a mutable reference to the owning Eigen::Matrix.
    MatrixType & container()
    {
      return m_matrix;
    }

    /// @brief Returns a const reference to the owning Eigen::Matrix.
    const MatrixType & container() const
    {
      return m_matrix;
    }

    /**
     * @brief Gets a const pointer to the beginning of the matrix's data buffer.
     *
     * @return A const pointer to the first element (`Scalar*`) of the matrix's raw data.
     *
     * @details This provides direct, low-level read-only access to the matrix's coefficients.
     *          The data is stored in a contiguous block, typically in column-major order for
     *          Eigen matrices. This is useful for interoperability with C-style APIs or other
     *          libraries that operate on raw memory buffers.
     *
     * @warning The returned pointer is only valid as long as this object exists and its
     *          `m_matrix` member is not reallocated. Accessing it after the object is
     *          destroyed leads to undefined behavior.
     */
    const Scalar * data() const
    {
      return m_matrix.data();
    }

    /**
     * @brief Gets a mutable pointer to the beginning of the map's data buffer.
     * @copydoc data() const
     * @return A mutable pointer to the first element (`Scalar*`) of the map's raw data.
     */
    Scalar * data()
    {
      return m_matrix.data();
    }

  }; // struct MatrixBlockElementTpl

  template<typename MatrixOrMap, typename Enable, typename DiagonalVectorType>
  BinaryOperator<
    internal::add_op,
    MatrixBlockElementTpl<MatrixOrMap, Enable>,
    Eigen::DiagonalWrapper<DiagonalVectorType>>
  operator+(
    const MatrixBlockElementTpl<MatrixOrMap, Enable> & matrix_block_elt,
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_matrix)
  {
    return {matrix_block_elt, diagonal_matrix};
  }

  template<typename MatrixOrMap, typename Enable, typename DiagonalVectorType>
  BinaryOperator<
    internal::sub_op,
    MatrixBlockElementTpl<MatrixOrMap, Enable>,
    Eigen::DiagonalWrapper<DiagonalVectorType>>
  operator-(
    const MatrixBlockElementTpl<MatrixOrMap, Enable> & matrix_block_elt,
    const Eigen::DiagonalWrapper<DiagonalVectorType> & diagonal_matrix)
  {
    return {matrix_block_elt, diagonal_matrix};
  }

} // namespace pinocchio
