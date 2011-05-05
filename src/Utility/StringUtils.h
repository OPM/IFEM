#pragma once
//==============================================================================
//!
//! \file StringUtils.h
//!
//! \date Feb 17 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Some general string manipulation functions
//!
//==============================================================================
#include <string>

std::string& replaceAll(std::string& context,
                        const std::string& from, const std::string& to);
