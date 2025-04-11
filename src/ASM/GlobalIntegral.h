// $Id$
//==============================================================================
//!
//! \file GlobalIntegral.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for classes representing integrated quantities.
//!
//==============================================================================

#ifndef _GLOBAL_INTEGRAL_H
#define _GLOBAL_INTEGRAL_H

#include <vector>
#include <cstddef>

class LocalIntegral;
struct Property;


/*!
  \brief Abstract base class representing a system level integrated quantity.
*/

class GlobalIntegral
{
public:
  //! \brief The default constructor.
  GlobalIntegral() {}
  //! \brief Empty destructor.
  virtual ~GlobalIntegral() {}

  //! \brief Initializes the integrated quantity to zero.
  virtual void initialize(char) {}
  //! \brief Finalizes the integrated quantity after element assembly.
  virtual bool finalize(bool) { return true; }

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  virtual bool assemble(const LocalIntegral*, int) { return true; }

  //! \brief Returns \e true if all elements can be assembled in parallel.
  virtual bool threadSafe() const { return false; }

  //! \brief Returns \e false if no contributions from a specified patch.
  virtual bool haveContributions(size_t, const std::vector<Property>&) const
  { return true; }
};

#endif
