// $Id$
//==============================================================================
//!
//! \file StringUtils.C
//!
//! \date Feb 17 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Some general string manipulation functions.
//!
//==============================================================================

#include "StringUtils.h"
#include <algorithm>


// taken from cppreference.com ::replace documentation
std::string& replaceAll(std::string& context, const std::string& from,
                        const std::string& to)
{
  size_t lookHere = 0;
  size_t foundHere = 0;
  while((foundHere = context.find(from, lookHere)) != std::string::npos)
  {
    context.replace(foundHere, from.size(), to);
    lookHere = foundHere + to.size();
  }
  return context;
}


std::vector<std::string> splitString(const std::string& str,
                                     int delimiter(int))
{
  std::vector<std::string> result;
  auto e=str.end();
  auto i=str.begin();
  while (i!=e) {
    i=std::find_if_not(i,e, delimiter);
    if (i==e)
      break;
    auto j=std::find_if(i,e, delimiter);
    result.push_back(std::string(i,j));
    i=j;
  }
  return result;
}
