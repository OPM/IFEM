// $Id$
//==============================================================================
//!
//! \file ASMutils.h
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utilities for assembly scope.
//!
//==============================================================================

#ifndef _ASM_UTILS_H
#define _ASM_UTILS_H

#include <iostream>
#include <vector>
#include <string>

class Vec3;

typedef std::vector<int>    IntVec; //!< General integer vector
typedef std::vector<IntVec> IntMat; //!< General 2D integer matrix


namespace ASM
{
  typedef std::pair<std::string,IntVec> NodeSet; //!< Named node set container

  //! \brief Creates a mesh by reading Matlab commands from an input stream.
  bool readMatlab(std::istream& is, IntMat& MNPC, std::vector<Vec3>& nodes,
                  std::vector<NodeSet>& nodeSets);
};

#endif
