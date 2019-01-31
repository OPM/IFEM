//==============================================================================
//!
//! \file readIntVec.h
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Read an integer array from a named file, for unit testing.
//!
//==============================================================================

#ifndef _READ_INT_VEC_H
#define _READ_INT_VEC_H

#include <fstream>
#include <vector>
#include <string>

typedef std::vector<int> IntVec; //!< General integer vector


/*!
  \brief Reads an integer array from a named file, for unit testing.
*/

static IntVec readIntVector (const std::string& file)
{
  std::ifstream f(file);
  size_t size;
  f >> size;
  if (!f) size = 0;
  IntVec result(size);
  for (int& i : result) f >> i;
  return result;
}

#endif
