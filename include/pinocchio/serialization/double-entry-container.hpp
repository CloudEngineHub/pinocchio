//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_serialization_double_entry_container_hpp__
#define __pinocchio_serialization_double_entry_container_hpp__

#include "pinocchio/serialization/eigen.hpp"
#include "pinocchio/container/double-entry-container.hpp"

namespace boost
{
  namespace serialization
  {

    namespace internal
    {
      template<typename StdVectorLike>
      struct DoubleEntryContainerAccessor
      : public ::pinocchio::container::DoubleEntryContainer<StdVectorLike>
      {
        typedef ::pinocchio::container::DoubleEntryContainer<StdVectorLike> Base;
        using Base::m_keys;
        using Base::m_values;
      };
    } // namespace internal

    template<typename Archive, typename StdVectorLike>
    void serialize(
      Archive & ar,
      ::pinocchio::container::DoubleEntryContainer<StdVectorLike> & container,
      const unsigned int /*version*/)
    {
      typedef internal::DoubleEntryContainerAccessor<StdVectorLike> Accessor;
      Accessor & container_ = static_cast<Accessor &>(container);
      ar & make_nvp("m_keys", container_.m_keys);
      ar & make_nvp("m_values", container_.m_values);
    }

  } // namespace serialization
} // namespace boost

#endif // ifndef __pinocchio_serialization_double_entry_container_hpp__
