//
// Copyright (c) 2015-2024 CNRS INRIA
//

#ifndef __pinocchio_collision_coal_convertion_hpp__
#define __pinocchio_collision_coal_convertion_hpp__

#include <coal/math/transform.h>
#include "pinocchio/spatial/se3.hpp"

namespace pinocchio
{
  template<typename Scalar>
  inline coal::Transform3s toCoalTransform3s(const SE3Tpl<Scalar> & m)
  {
    SE3Tpl<double, 0> m_ = m.template cast<double>();
    return coal::Transform3s(m_.rotation(), m_.translation());
  }

  inline SE3 toPinocchioSE3(const coal::Transform3s & tf)
  {
    typedef SE3::Scalar Scalar;
    return SE3(tf.getRotation().cast<Scalar>(), tf.getTranslation().cast<Scalar>());
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_collision_coal_convertion_hpp__
