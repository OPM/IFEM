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
#include <iostream>
#include <vector>


/*!
  \brief Base class for NURBS-based FEM simulators with input file parsing.
*/

class SIMinput
{
protected:
  //! \brief The constructor initializes \a myPid, \a nProc and \a myHeading.
  SIMinput(const char* heading = NULL);
  //! \brief Constructor using control parameters from another SIMinput object.
  SIMinput(SIMinput& anotherSIM);

public:
  //! \brief Empty destructor.
  virtual ~SIMinput() {}

  //! \brief Reads model data from the specified input file \a *fileName.
  virtual bool read(const char* fileName);

private:
  //! \brief Reads a flat text input file.
  bool readFlat(const char* fileName);
  //! \brief Reads an XML input file.
  bool readXML(const char* fileName);

protected:
  //! \brief Handles the parsing order for certain XML-tags.
  //! \param[in] base The base tag containing the elements to be prioritized
  //! \param[out] parsed Vector of XML-elements that was parsed
  //!
  //! \details Certain tags need to be parsed before others. This method takes
  //! care of this. It is called by the \a readXML method in order to read the
  //! top level tags in the required order. It can also be called by the
  //! application-specific SIM class prior to parsing its data blocks.
  //! In that case the \a getPrioritizedTags method should be reimplemented
  //! by the sub-class to take care of the application-specific tags.
  bool handlePriorityTags(const TiXmlElement* base,
			  std::vector<const TiXmlElement*>& parsed);

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const { return NULL; }

public:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is) = 0;
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
  std::string inputFile; //!< The parsed input file
};

#endif
