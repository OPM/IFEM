//==============================================================================
//!
//! \file StringUtils.C
//!
//! \date Feb 17 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Some general string manipulation functions
//!
//==============================================================================

#include <StringUtils.h>

// taken from cppreference.com ::replace documentation
std::string& replaceAll(std::string& context, const std::string& from,
                        const std::string& to)
{
  size_t lookHere = 0;
  size_t foundHere;
  while((foundHere = context.find(from, lookHere)) != std::string::npos)
  {
    context.replace(foundHere, from.size(), to);
    lookHere = foundHere + to.size();
  }
  return context;
}
