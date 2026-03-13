//
// Copyright (c) 2026 INRIA
//
#pragma once

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/serialization.hpp"
#endif // PINOCCHIO_LSP

#define BOOST_SERIALIZATION_MAKE_NVP(member) boost::serialization::make_nvp(##member, member)
