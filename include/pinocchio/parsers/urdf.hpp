//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <cassert>
#include <iosfwd>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <stdexcept>

#include <boost/optional.hpp>
#include <boost/none.hpp>

#include <Eigen/Core>

#include <urdf_model/model.h>
#include <urdf_model/pose.h>

#include "pinocchio/macros.hpp"

#include "pinocchio/spatial.hpp"
#include "pinocchio/multibody.hpp"

#include "pinocchio/parsers/config.hpp"
#include "pinocchio/parsers/fwd.hxx"
#include "pinocchio/parsers/meshloader-fwd.hxx"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/parsers/urdf/urdf.hxx"
#include "pinocchio/parsers/urdf/utils.hxx"
#include "pinocchio/parsers/urdf/model.hxx"
#include "pinocchio/parsers/urdf/geometry.hxx"
// IWYU pragma: end_exports
