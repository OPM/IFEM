// $Id$
//==============================================================================
//!
//! \file SIMargsBase.h
//!
//! \date Jul 15 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for pre-parsing of XML input files for simulators.
//!
//==============================================================================

#ifndef _SIM_ARGS_BASE_H_
#define _SIM_ARGS_BASE_H_

#include "XMLInputBase.h"


/*!
  \brief Base class for input file pre-parsing in applications.
*/

class SIMargsBase : public XMLInputBase
{
public:
  //! \brief The constructor initializes the default parameter values.
  explicit SIMargsBase(const char* ctx) : context(ctx), dim(3), adap(false) {}

  //! \brief Parses a command-line argument.
  virtual bool parseArg(const char* argv);

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

private:
  const char* context; //!< Application-specific context tag

public:
  int  dim;  //!< Dimensionality of simulation
  bool adap; //!< If \e true, run an adaptive simulator
};

#endif
