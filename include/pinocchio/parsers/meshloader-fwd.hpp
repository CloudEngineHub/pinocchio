//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_parsers_meshloader_fwd_hpp__
#define __pinocchio_parsers_meshloader_fwd_hpp__

#include <memory>

#ifdef PINOCCHIO_WITH_COLLISION
  #include <coal/config.hh>
#endif // PINOCCHIO_WITH_COLLISION

namespace coal
{
  class MeshLoader;
  typedef std::shared_ptr<MeshLoader> MeshLoaderPtr;
} // namespace coal

#endif // __pinocchio_parsers_meshloader_fwd_hpp__
