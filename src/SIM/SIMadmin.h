// $Id$
//==============================================================================
//!
//! \file SIMadmin.h
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration base class for FEM simulators.
//!
//==============================================================================

#ifndef _SIM_ADMIN_H
#define _SIM_ADMIN_H

#include "XMLInputBase.h"
#include "SIMoptions.h"
#include "ProcessAdm.h"
#include <iostream>
#include <string>


/*!
  \brief Administration base class for FEM simulators.
  \details This class serves as a common base for all types of simulator drivers
  in IFEM and contains the general top-level methods for reading the model input
  file, as well as data for administration of parallel executions.
*/

class SIMadmin : public XMLInputBase
{
protected:
  //! \brief The default constructor initializes the process administrator.
  explicit SIMadmin(const char* heading = nullptr);
  //! \brief Copy constructor.
  SIMadmin(SIMadmin& anotherSIM);

public:
  //! \brief Empty destructor.
  virtual ~SIMadmin() {}

  //! \brief Reads model data from the specified input file \a *fileName.
  virtual bool read(const char* fileName);

  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual bool preprocess(const std::vector<int>&, bool) { return false; }

  //! \brief Returns the parallel process administrator.
  const ProcessAdm& getProcessAdm() const { return adm; }

  //! \brief Returns the global process ID.
  //! \note May be different from the process ID used in the equation solver.
  int getGlobalProcessID() const { return myPid; }

  //! \brief Returns the simulator heading.
  const std::string& getHeading() const { return myHeading; }
  //! \brief Defines the simulator heading.
  void setHeading(const std::string& heading) { myHeading = heading; }
  //! \brief Prints the heading of this simulator, if any, to IFEM::cout.
  void printHeading(int& supStep) const;

private:
  //! \brief Reads a flat text input file (the old file format).
  bool readFlat(const char* fileName);

public:
  static int  msgLevel;  //!< Controls the console output amount during solving
  SIMoptions& opt;       //!< Simulation control parameters

private:
  SIMoptions  myOpts;    //!< Actual control parameters owned by this simulator

protected:
  ProcessAdm  adm;       //!< Parallel administrator
  int         myPid;     //!< Processor ID in parallel simulations
  int         nProc;     //!< Number of processors in parallel simulations
  std::string myHeading; //!< Heading written before reading the input file
};

#endif
