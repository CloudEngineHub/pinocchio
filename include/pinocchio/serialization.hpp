//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <new>
#include <cstddef>
#include <vector>
#include <fstream>
#include <locale>
#include <string>
#include <sstream>
#include <stdexcept>
#include <array>

#include <boost/version.hpp>
#if BOOST_VERSION / 100 % 1000 == 78 && __APPLE__
// See https://github.com/qcscine/utilities/issues/5#issuecomment-1246897049 for further details

  #ifndef BOOST_ASIO_DISABLE_STD_ALIGNED_ALLOC
    #define DEFINE_BOOST_ASIO_DISABLE_STD_ALIGNED_ALLOC
    #define BOOST_ASIO_DISABLE_STD_ALIGNED_ALLOC
  #endif

  #include <boost/asio/streambuf.hpp>

  #ifdef DEFINE_BOOST_ASIO_DISABLE_STD_ALIGNED_ALLOC
    #undef BOOST_ASIO_DISABLE_STD_ALIGNED_ALLOC
  #endif

#else
  #include <boost/asio/streambuf.hpp>
#endif
#include <boost/archive/basic_archive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/throw_exception.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/blank.hpp>

#include <Eigen/Core>

#ifdef PINOCCHIO_WITH_COLLISION
  #define COAL_SKIP_EIGEN_BOOST_SERIALIZATION
  #include <coal/serialization/collision_data.h>
  #undef COAL_SKIP_EIGEN_BOOST_SERIALIZATION
  #include <coal/serialization/geometric_shapes.h>
  #include <coal/serialization/hfield.h>
  #include <coal/serialization/octree.h>
  #include <coal/serialization/convex.h>
  #include <coal/serialization/BVH_model.h>
#endif // PINOCCHIO_WITH_COLLISION

#include "pinocchio/macros.hpp"

#include "pinocchio/container/double-entry-container.hpp"
#include "pinocchio/container/matrix-stack.hpp"
#include "pinocchio/container/eigen-storage.hxx"

#include "pinocchio/math.hpp"
#include "pinocchio/spatial.hpp"
#include "pinocchio/multibody.hpp"
#include "pinocchio/multibody/joint.hpp"
#include "pinocchio/geometry.hpp"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/serialization/fwd.hxx"
#include "pinocchio/serialization/common.hxx"
#include "pinocchio/serialization/static-buffer.hxx"
#include "pinocchio/serialization/boost-blank.hxx"
#include "pinocchio/serialization/archive.hxx"
#include "pinocchio/serialization/serializable.hxx"
#include "pinocchio/serialization/vector.hxx"
#include "pinocchio/serialization/coal.hxx"
#include "pinocchio/serialization/csv.hxx"
#include "pinocchio/serialization/eigen.hxx"
#include "pinocchio/serialization/matrix-stack.hxx"
#include "pinocchio/serialization/double-entry-container.hxx"
#include "pinocchio/serialization/eigen-storage.hxx"
#include "pinocchio/serialization/symmetric3.hxx"
#include "pinocchio/serialization/se3.hxx"
#include "pinocchio/serialization/motion.hxx"
#include "pinocchio/serialization/force.hxx"
#include "pinocchio/serialization/inertia.hxx"
#include "pinocchio/serialization/joints-motion-subspace.hxx"
#include "pinocchio/serialization/joints-motion.hxx"
#include "pinocchio/serialization/joints-transform.hxx"
#include "pinocchio/serialization/joints-data.hxx"
#include "pinocchio/serialization/joints-model.hxx"
#include "pinocchio/serialization/frame.hxx"
#include "pinocchio/serialization/model.hxx"
#include "pinocchio/serialization/data.hxx"
// IWYU pragma: end_exports
