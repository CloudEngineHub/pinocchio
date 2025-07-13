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
    typedef EigenStorageTpl<NewPlainMatrixType> type;
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

    /// \brief Default constructor
    ///
    /// \param[in] max_elts Maximum number of matrices contained in the stack
    /// \param[in] max_elt_size Maximal size of each matrices (rows() x cols()) if known at
    /// construction time. Default value to 0.
    ///
    MatrixStackTpl(const std::size_t max_elts, const std::size_t max_elt_size = 0)
    : offsets()
    , data_ptr(nullptr)
    , memory_capacity(0)
    {
      if (max_elts > 0)
      {
        offsets.reserve(max_elts);
        matrix_maps.reserve(max_elts);

        // Allocate the full memory if max_elt_size is given
        if (max_elt_size > 0)
        {
          const std::size_t max_chunck_size = max_elt_size * sizeof(Scalar) + Alignment;
          const std::size_t max_total_size = max_elts * max_chunck_size;

          data_ptr =
            MatrixStackTpl::malloc(max_total_size); // the first element is for sure aligned
          memory_capacity = data_ptr != nullptr ? max_total_size : 0;
        }
      }
    }

    void push_back(const Index rows, const Index cols)
    {
      void * next_data_ptr = matrix_maps.size() == 0
                               ? data_ptr
                               : inc_ptr(matrix_maps.back().data(), raw_size(matrix_maps.back()));
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

      const std::size_t matrix_raw_size = std::size_t(rows * cols) * sizeof(Scalar);
      const std::size_t loss_bits =
        (reinterpret_cast<std::size_t>(aligned_data)
         - reinterpret_cast<std::size_t>(next_data_ptr));
      const std::size_t new_memory_chunck_size = matrix_raw_size + loss_bits;

      const std::size_t current_memory_size =
        reinterpret_cast<std::size_t>(next_data_ptr) - reinterpret_cast<std::size_t>(data_ptr);
      if (current_memory_size + new_memory_chunck_size > memory_capacity)
      { // We need to proceed to a new allocation
        const std::size_t new_size =
          2 * (current_memory_size + new_memory_chunck_size); // we double the allocated chunck

        if (data_ptr == nullptr)
        {
          data_ptr = MatrixStackTpl::malloc(new_size);
        }
        else
        {
          data_ptr = MatrixStackTpl::realloc(data_ptr, new_size, memory_capacity);
        }
        assert(data_ptr != nullptr);
        memory_capacity = new_size;

        // We need to realign all the existing Eigen maps
        for (std::size_t i = 0; i < matrix_maps.size(); ++i)
        {
          auto & matrix_map = matrix_maps[i];
          const auto offset_value = offsets[i];

          void * new_map_data_ptr = inc_ptr(data_ptr, offset_value);

          new (&matrix_map) MapType(
            reinterpret_cast<Scalar *>(new_map_data_ptr), matrix_map.rows(), matrix_map.cols());
        }

        void * next_data_ptr = matrix_maps.size() == 0
                                 ? data_ptr
                                 : inc_ptr(matrix_maps.back().data(), raw_size(matrix_maps.back()));
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
      matrix_maps.push_back(aligned_map);

      offsets.push_back(
        reinterpret_cast<std::size_t>(aligned_data) - reinterpret_cast<std::size_t>(data_ptr));
    }

    /// \brief Returns a reference to the last element in the container.
    RefMapType back()
    {
      return matrix_maps.back();
    }
    /// \brief Returns a reference to the last element in the container.
    ConstRefMapType back() const
    {
      return matrix_maps.back();
    }

    ///  \brief Checks if the container has no elements.
    ///
    ///  \returns true if the container is empty, false otherwise.
    bool empty() const
    {
      return matrix_maps.empty();
    }

    std::size_t capacity() const
    {
      return matrix_maps.capacity();
    }

    /// \brief Returns a reference to the element at specified location pos.
    RefMapType operator[](const std::size_t pos)
    {
      return matrix_maps[pos];
    }
    /// \brief Returns a reference to the element at specified location pos.
    ConstRefMapType operator[](const std::size_t pos) const
    {
      return matrix_maps[pos];
    }

    /// \brief Returns the number of elements in the container.
    std::size_t size() const
    {
      return matrix_maps.size();
    }

    /// \brief Returns a pointer to the underlying array serving as element storage.
    void * data()
    {
      return data_ptr;
    }
    /// \brief Returns a pointer to the underlying array serving as element storage.
    const void * data() const
    {
      return data_ptr;
    }

    ~MatrixStackTpl()
    {
      MatrixStackTpl::free(data_ptr);
    }

  protected:
    static void * malloc(std::size_t size, std::size_t alignment = Alignment)
    {
      return Eigen::internal::handmade_aligned_malloc(
        size - alignment, /* we can remove one alignment value already taken into account in
                             max_total_size*/
        alignment);
    }

    static void free(void * ptr)
    {
      Eigen::internal::handmade_aligned_free(ptr);
    }

    static void * realloc(
      void * ptr, std::size_t new_size, std::size_t old_size, std::size_t alignment = Alignment)
    {
      return Eigen::internal::handmade_aligned_realloc(ptr, new_size, old_size, alignment);
    }

    static void * inc_ptr(void * ptr, std::size_t inc_value)
    {
      return reinterpret_cast<void *>(reinterpret_cast<std::size_t>(ptr) + inc_value);
    }

    static std::size_t raw_size(const MapType & map)
    {
      return sizeof(Scalar) * std::size_t(map.rows() * map.cols());
    }

    std::vector<std::size_t> offsets;
    std::vector<MapType> matrix_maps;
    void * data_ptr;
    std::size_t memory_capacity;
  }; // struct MatrixStackTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_container_matrix_stack_hpp__
