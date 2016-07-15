// $Id$
//==============================================================================
//!
//! \file SIMinput.h
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for simulators with input parsing functionality.
//!
//==============================================================================

#ifndef _SIM_INPUT_H
#define _SIM_INPUT_H

#include "SIMoptions.h"
#include "ProcessAdm.h"
#include "XMLInputBase.h"
#include <iostream>
#include <vector>


/*!
  \brief Base class for NURBS-based FEM simulators with input file parsing.
*/

class SIMinput : public XMLInputBase
{
protected:
  //! \brief The default constructor initializes the process administrator.
  SIMinput(const char* heading = nullptr);
  //! \brief Copy constructor.
  SIMinput(SIMinput& anotherSIM);

public:
  //! \brief Empty destructor.
  virtual ~SIMinput() {}

  //! \brief Reads model data from the specified input file \a *fileName.
  virtual bool read(const char* fileName);

public:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) = 0;

  //! \brief Returns the parallel process administrator.
  const ProcessAdm& getProcessAdm() const { return adm; }

  //! \brief Returns the global process ID.
  //! \note May be different from the process ID used in the equation solver.
  int getGlobalProcessID() const { return myPid; }

  //! \brief Returns the simulator heading.
  const std::string& getHeading() const { return myHeading; }
  //! \brief Defines the simulator heading.
  void setHeading(const std::string& heading) { myHeading = heading; }

protected:
  //! \brief Prints the heading of this (sub-step) solver, if any, to std::cout.
  void printHeading(int& supStep) const;

  //! \brief Reads a flat text input file (obsolete file format).
  bool readFlat(const char* fileName);

public:
  static int msgLevel; //!< Controls the amount of console output during solving
  SIMoptions& opt;     //!< Simulation control parameters

private:
  SIMoptions myOpts; //!< The actual control parameters owned by this simulator

protected:
  ProcessAdm adm; //!< Parallel administrator
  int myPid;      //!< Processor ID in parallel simulations
  int nProc;      //!< Number of processors in parallel simulations

  std::string myHeading; //!< Heading written before reading the input file
};

#endif
