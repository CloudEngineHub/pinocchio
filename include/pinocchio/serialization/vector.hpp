//
// Copyright (c) 2019-2021 INRIA
//

#ifndef __pinocchio_serialization_vector_hpp__
#define __pinocchio_serialization_vector_hpp__

#include <vector>

#include <boost/version.hpp>
#include <boost/core/addressof.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

namespace boost
{
  namespace serialization
  {

    template<class T, class Allocator>
    inline const nvp<std::vector<T, Allocator>>
    make_nvp(const char * name, std::vector<T, Allocator> & t)
    {
      return nvp<std::vector<T, Allocator>>(name, t);
    }

  } // namespace serialization
} // namespace boost

#endif // ifndef __pinocchio_serialization_vector_hpp__
