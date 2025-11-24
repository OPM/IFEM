// $Id$
//==============================================================================
//!
//! \file Catch2Support.h
//!
//! \date Sep 9 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief IFEM catch2 support
//!
//==============================================================================

#ifndef CATCH2_SUPPORT_H_
#define CATCH2_SUPPORT_H_

#if CATCH2_VERSION_MAJOR > 2
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#else
#include <catch2/catch.hpp>
#endif

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

#endif
