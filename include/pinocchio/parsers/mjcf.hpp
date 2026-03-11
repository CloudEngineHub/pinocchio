//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep

#include <cstddef>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <iostream>

#include <Eigen/Core>

#include <boost/math/constants/constants.hpp>
#include <boost/none.hpp>
#include <boost/none_t.hpp>
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant.hpp>
#include <boost/logic/tribool.hpp>

#include "pinocchio/spatial.hpp"

#include "pinocchio/multibody.hpp"

#include "pinocchio/geometry.hpp"

#include "pinocchio/algorithm/constraints.hpp"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/meshloader-fwd.hxx"
#include "pinocchio/parsers/config.hpp"
#include "pinocchio/parsers/scalar-model.hxx"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/parsers/mjcf/mjcf.hxx"
#include "pinocchio/parsers/mjcf/mjcf-graph.hxx"
#include "pinocchio/parsers/mjcf/model.hxx"
#include "pinocchio/parsers/mjcf/geometry.hxx"
// IWYU pragma: end_exports
