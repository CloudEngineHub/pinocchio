//
// Copyright (c) 2026 INRIA
//

#ifndef __pinocchio_collision_fcl_convertion_hpp__
#define __pinocchio_collision_fcl_convertion_hpp__

#include "pinocchio/macros.hpp"

// clang-format off
PINOCCHIO_MOVED_HEADER_PINOCCHIO4(pinocchio/collision/fcl-pinocchio-conversions.hpp, pinocchio/collision/coal-pinocchio-conversions.hpp)
// clang-format on

#include "pinocchio/collision/coal-pinocchio-conversions.hpp"

namespace pinocchio
{

  template<typename Scalar>
  PINOCCHIO_DEPRECATED inline coal::Transform3s toFclTransform3f(const SE3Tpl<Scalar> & m)
  {
    return toCoalTransform3s(m);
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_collision_fcl_convertion_hpp__
