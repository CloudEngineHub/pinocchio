//
// Copyright (c) 2026 INRIA
//
#pragma once

// IWYU pragma: begin_keep
#ifdef WIN32
  #include <Windows.h>
  #include <stdint.h> // portable: uint64_t   MSVC: __int64
#else
  #include <sys/time.h>
#endif
#include <iostream>
#include <stack>
// IWYU pragma: end_keep

// IWYU pragma: begin_exports
#include "pinocchio/src/utils/timer.hxx"
// IWYU pragma: export
