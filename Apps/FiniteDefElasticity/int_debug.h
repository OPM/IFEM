// $Id$
//==============================================================================
//!
//! \file int_debug.h
//!
//! \date Jul 10 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Convenience define for debugging of finite deformation integrands.
//!
//==============================================================================

#ifdef USE_OPENMP
#if INT_DEBUG > 0
#undef INT_DEBUG
#endif
#endif

#ifdef USE_FTNMAT
#ifndef INT_DEBUG
#define INT_DEBUG 0
#endif
#endif
