// $Id$
//==============================================================================
//!
//! \file StringUtils.h
//!
//! \date Feb 17 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Some general string manipulation functions.
//!
//==============================================================================

#pragma once
#include <string>
#include <vector>

//! \brief Replaces all occurances of \a from in \a context to \a to.
std::string& replaceAll(std::string& context,
                        const std::string& from, const std::string& to);

//! \brief Split a string on a given delimiter
//! \param[in] str The string to split
//! \param[in] delimiter Function callback for the delimiter to use
std::vector<std::string> splitString(const std::string& str,
                                     int delimiter(int) = ::isspace);
