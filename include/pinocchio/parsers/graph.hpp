//
// Copyright (c) INRIA 2026
//
#pragma once

// IWYU pragma: begin_keep
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <utility>

#include <boost/optional.hpp>
#include <boost/none.hpp>
#include <boost/variant.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "pinocchio/context.hpp"
#include "pinocchio/macros.hpp"

#include "pinocchio/utils/check.hpp"

#include "pinocchio/spatial.hpp"
#include "pinocchio/math.hpp"
#include "pinocchio/multibody.hpp"
#include "pinocchio/multibody/joint.hpp"

#include "pinocchio/parsers/config.hpp"
#include "pinocchio/parsers/meshloader-fwd.hxx"
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/parsers/graph/fwd.hxx"
#include "pinocchio/parsers/graph/frames.hxx"
#include "pinocchio/parsers/graph/joints.hxx"
#include "pinocchio/parsers/graph/geometries.hxx"
#include "pinocchio/parsers/graph/model-graph.hxx"
#include "pinocchio/parsers/graph/graph-visitor.hxx"
#include "pinocchio/parsers/graph/model-graph-algo.hxx"
#include "pinocchio/parsers/graph/model-graph-algo-geometry.hxx"
#include "pinocchio/parsers/graph/model-configuration-converter.hxx"
// IWYU pragma: end_exports
