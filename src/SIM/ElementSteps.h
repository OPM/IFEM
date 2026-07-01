// $Id$
//==============================================================================
//!
//! \file ElementSteps.h
//!
//! \date Jun 30 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Spatial element step function.
//!
//==============================================================================

#ifndef _ELEMENT_STEPS_H
#define _ELEMENT_STEPS_H

#include "FunctionSum.h"

class SIMbase;


/*!
  \brief A sum of single-element step functions.
*/

class ElementSteps : public RealFuncSum
{
public:
  //! \brief The constructor reads the function definition from a string.
  //! \param[in] input The string to parse for function parameters
  //! \param[in] sim The simulator holding the element information
  //! \param[in] nsd Number of spatial dimensions
  ElementSteps(const char* input, const SIMbase& sim, int nsd);
};

#endif
