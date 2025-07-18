//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_container_matrix_stack_hpp__
#define __pinocchio_container_matrix_stack_hpp__

#include "pinocchio/fwd.hpp"

namespace pinocchio
{

  template<typename MatrixLike, std::size_t Alignment = EIGEN_DEFAULT_ALIGN_BYTES>
  struct MatrixStackTpl;

  template<typename NewScalar, typename MatrixLike, std::size_t Alignment>
  struct CastType<NewScalar, MatrixStackTpl<MatrixLike, Alignment>>
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) PlainMatrixType;

    typedef typename PlainMatrixType::template CastXpr<NewScalar>::Type NewPlainMatrixExpression;
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(
      typename std::remove_reference<
        typename std::remove_const<NewPlainMatrixExpression>::type>::type) NewPlainMatrixType;
    typedef MatrixStackTpl<NewPlainMatrixType, Alignment> type;
  };

  template<typename MatrixLike, std::size_t _Alignment>
  struct MatrixStackTpl
  {
    typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) PlainMatrixType;
    typedef typename PlainMatrixType::Scalar Scalar;
    typedef typename Eigen::Index Index;

    enum
    {
      Alignment = _Alignment
    };

    typedef Eigen::Map<PlainMatrixType, Alignment> MapType;
    typedef MapType & RefMapType;
    typedef const MapType & ConstRefMapType;

    typedef const Eigen::Map<const PlainMatrixType, Alignment> ConstMapType;
    typedef ConstMapType & RefConstMapType;

    typedef std::vector<MapType> MapVector;

    typedef typename MapVector::iterator iterator;
    typedef typename MapVector::const_iterator const_iterator;

    /// \brief Default constructor
    MatrixStackTpl()
    : m_data_ptr(nullptr)
    , m_memory_capacity(0)
    {
    }

    /// \brief Constructor
    ///
    /// \param[in] max_elts Maximum number of matrices contained in the stack
    /// \param[in] max_elt_size Maximal size of each matrices (rows() x cols()) if known at
    /// construction time. Default value to 0.
    ///
    explicit MatrixStackTpl(const std::size_t max_elts, const std::size_t max_elt_size = 0)
    : m_offsets()
    , m_data_ptr(nullptr)
    , m_memory_capacity(0)
    {
      if (max_elts > 0)
      {
        m_offsets.reserve(max_elts);
        m_matrix_maps.reserve(max_elts);

        // Allocate the full memory if max_elt_size is given
        if (max_elt_size > 0)
        {
          const std::size_t max_chunck_size = max_elt_size * sizeof(Scalar);
          const std::size_t max_total_size =
            max_elts * max_chunck_size + (max_elt_size - 1) * Alignment;

          m_data_ptr =
            MatrixStackTpl::malloc(max_total_size); // the first element is for sure aligned
          m_memory_capacity = m_data_ptr != nullptr ? max_total_size : 0;
        }
      }
    }

    ///
    /// @brief Copy constructor. Creates a deep copy of *this.
    /// @param other MatrixStackTpl to copy
    ///
    MatrixStackTpl(const MatrixStackTpl & other)
    : m_data_ptr(nullptr)
    {
      *this = other;
    }

    ///
    /// @brief Copy operator
    /// @param other MatrixStackTpl to copy
    /// @returns a reference to this.
    ///
    MatrixStackTpl & operator=(const MatrixStackTpl & other)
    {
      free(m_data_ptr);

      m_memory_capacity = other.raw_size();

      if (m_memory_capacity == 0)
      {
        m_data_ptr = nullptr;
        return *this;
      }

      m_data_ptr = MatrixStackTpl::malloc(m_memory_capacity);
      if (m_data_ptr == nullptr)
      {
        m_memory_capacity = 0;
        return *this;
      }

      // Copy raw data
      std::memcpy(m_data_ptr, other.m_data_ptr, m_memory_capacity);

      // Add aligned map
      m_matrix_maps.clear();
      m_matrix_maps.reserve(other.m_matrix_maps.size());
      m_offsets = other.m_offsets;
      for (std::size_t i = 0; i < other.m_matrix_maps.size(); ++i)
      {
        const auto offset_value = m_offsets[i];
        const auto & other_matrix_map = other.m_matrix_maps[i];

        void * aligned_data = incr_ptr(m_data_ptr, offset_value);
        assert(
          reinterpret_cast<std::size_t>(aligned_data) % Alignment == 0
          && "aligned_data is not properly aligned.");

        MapType aligned_map = MapType(
          reinterpret_cast<Scalar *>(aligned_data), other_matrix_map.rows(),
          other_matrix_map.cols());
        // aligned_map = other_matrix_map; // copy data
        m_matrix_maps.push_back(aligned_map);
      }

      return *this;
    }

    /// @brief Equality comparison operator.
    /// @param other MatrixStackTpl to compare with.
    /// @returns true if the underlying maps are equal.
    bool operator==(const MatrixStackTpl & other) const
    {
      if (this == &other)
        return true;
      return m_matrix_maps == other.m_matrix_maps;
    }

    /// @brief Inequality comparison operator.
    /// @param other MatrixStackTpl to compare with.
    /// @return true if the underlying maps are not equal.
    bool operator!=(const MatrixStackTpl & other) const
    {
      return !(*this == other);
    }

    void push_back(const Index rows, const Index cols)
    {
      void * next_data_ptr =
        m_matrix_maps.size() == 0
          ? m_data_ptr
          : incr_ptr(m_matrix_maps.back().data(), raw_map_size(m_matrix_maps.back()));
      void * aligned_data =
        reinterpret_cast<std::size_t>(next_data_ptr) % Alignment == 0
          ? /* next_data_ptr is aligned */
          next_data_ptr
          : reinterpret_cast<void *>(
              (reinterpret_cast<std::size_t>(next_data_ptr) & ~(std::size_t(Alignment - 1)))
              + Alignment);
      assert(
        reinterpret_cast<std::size_t>(aligned_data) % Alignment == 0
        && "aligned_data is not properly aligned.");

      const std::size_t matrix_raw_map_size = std::size_t(rows * cols) * sizeof(Scalar);
      const std::size_t loss_bits =
        (reinterpret_cast<std::size_t>(aligned_data)
         - reinterpret_cast<std::size_t>(next_data_ptr));
      const std::size_t new_memory_chunck_size = matrix_raw_map_size + loss_bits;

      const std::size_t current_memory_size =
        reinterpret_cast<std::size_t>(next_data_ptr) - reinterpret_cast<std::size_t>(m_data_ptr);
      if (current_memory_size + new_memory_chunck_size > m_memory_capacity)
      { // We need to proceed to a new allocation
        const std::size_t new_size =
          2 * (current_memory_size + new_memory_chunck_size); // we double the allocated chunck

        if (m_data_ptr == nullptr)
        {
          m_data_ptr = MatrixStackTpl::malloc(new_size);
        }
        else
        {
          m_data_ptr = MatrixStackTpl::realloc(m_data_ptr, new_size, m_memory_capacity);
        }
        assert(m_data_ptr != nullptr);
        m_memory_capacity = new_size;

        // We need to realign all the existing Eigen maps
        for (std::size_t i = 0; i < m_matrix_maps.size(); ++i)
        {
          auto & matrix_map = m_matrix_maps[i];
          const auto offset_value = m_offsets[i];

          void * new_map_data_ptr = incr_ptr(m_data_ptr, offset_value);

          new (&matrix_map) MapType(
            reinterpret_cast<Scalar *>(new_map_data_ptr), matrix_map.rows(), matrix_map.cols());
        }

        void * next_data_ptr =
          m_matrix_maps.size() == 0
            ? m_data_ptr
            : incr_ptr(m_matrix_maps.back().data(), raw_map_size(m_matrix_maps.back()));
        aligned_data =
          reinterpret_cast<std::size_t>(next_data_ptr) % Alignment == 0
            ? /* next_data_ptr is aligned */
            next_data_ptr
            : reinterpret_cast<void *>(
                (reinterpret_cast<std::size_t>(next_data_ptr) & ~(std::size_t(Alignment - 1)))
                + Alignment);
        assert(
          reinterpret_cast<std::size_t>(aligned_data) % Alignment == 0
          && "aligned_data is not properly aligned.");
      }

      MapType aligned_map = MapType(reinterpret_cast<Scalar *>(aligned_data), rows, cols);
      m_matrix_maps.push_back(aligned_map);

      m_offsets.push_back(
        reinterpret_cast<std::size_t>(aligned_data) - reinterpret_cast<std::size_t>(m_data_ptr));
    }

    /// \brief Returns a reference to the last element in the container.
    RefMapType back()
    {
      return m_matrix_maps.back();
    }
    /// \brief Returns a reference to the last element in the container.
    ConstRefMapType back() const
    {
      return m_matrix_maps.back();
    }

    ///  \brief Checks if the container has no elements.
    ///
    ///  \returns true if the container is empty, false otherwise.
    bool empty() const
    {
      return m_matrix_maps.empty();
    }

    std::size_t capacity() const
    {
      return m_matrix_maps.capacity();
    }

    /// \brief Returns a reference to the element at specified location pos.
    RefMapType operator[](const std::size_t pos)
    {
      return m_matrix_maps[pos];
    }
    /// \brief Returns a reference to the element at specified location pos.
    ConstRefMapType operator[](const std::size_t pos) const
    {
      return m_matrix_maps[pos];
    }

    /// \brief Returns the number of elements in the container.
    std::size_t size() const
    {
      return m_matrix_maps.size();
    }

    /// \brief Returns a pointer to the underlying array serving as element storage.
    void * data()
    {
      return m_data_ptr;
    }
    /// \brief Returns a pointer to the underlying array serving as element storage.
    const void * data() const
    {
      return m_data_ptr;
    }

    void clear()
    {
      m_offsets.clear();
      m_matrix_maps.clear();
    }

    iterator begin()
    {
      return m_matrix_maps.begin();
    }

    iterator end()
    {
      return m_matrix_maps.end();
    }

    const_iterator begin() const
    {
      return m_matrix_maps.begin();
    }

    const_iterator end() const
    {
      return m_matrix_maps.end();
    }

    iterator rbegin()
    {
      return m_matrix_maps.cbegin();
    }

    iterator rend()
    {
      return m_matrix_maps.cend();
    }

    const_iterator rbegin() const
    {
      return m_matrix_maps.cbegin();
    }

    const_iterator rend() const
    {
      return m_matrix_maps.cend();
    }

    void apply(const std::function<void(MapType)> & func)
    {
      std::for_each(begin(), end(), func);
    }

    void apply(const std::function<void(const MapType)> & func) const
    {
      std::for_each(begin(), end(), func);
    }

    ~MatrixStackTpl()
    {
      MatrixStackTpl::free(m_data_ptr);
    }

  protected:
    static void * malloc(std::size_t size, std::size_t alignment = Alignment)
    {
      assert(size > 0 && "size should be greater than 0.");
      return Eigen::internal::handmade_aligned_malloc(size, alignment);
    }

    static void free(void * ptr)
    {
      Eigen::internal::handmade_aligned_free(ptr);
    }

    static void * realloc(
      void * ptr, std::size_t new_size, std::size_t old_size, std::size_t alignment = Alignment)
    {
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
      return Eigen::internal::handmade_aligned_realloc(ptr, new_size, old_size, alignment);
#else
      return Eigen::internal::handmade_aligned_realloc(ptr, new_size, old_size);
#endif
    }

    static void * incr_ptr(void * ptr, std::size_t inc_value)
    {
      return reinterpret_cast<void *>(reinterpret_cast<std::size_t>(ptr) + inc_value);
    }

    static std::size_t raw_map_size(const MapType & map)
    {
      return sizeof(Scalar) * std::size_t(map.rows() * map.cols());
    }

    std::size_t raw_size() const
    {
      return m_matrix_maps.size() == 0
               ? 0
               : m_offsets.back() + raw_map_size(m_matrix_maps.back()); // + Alignment;
    }

    std::vector<std::size_t> m_offsets;
    MapVector m_matrix_maps;
    void * m_data_ptr;
    std::size_t m_memory_capacity;
  }; // struct MatrixStackTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_container_matrix_stack_hpp__
