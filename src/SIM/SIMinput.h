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

#include <iostream>
#include <vector>


/*!
  \brief Base class for NURBS-based FEM simulators with input file parsing.
*/

class TiXmlElement;

class SIMinput
{
protected:
  //! \brief The default constructor initializes \a myPid and \a nProc.
  SIMinput();

public:
  //! \brief Empty destructor.
  virtual ~SIMinput() {}

  //! \brief Reads model data from the specified input file \a *fileName.
  virtual bool read(const char* fileName);
  //! \brief Reads a flat text input file
  virtual bool readFlat(const char* fileName);
  //! \brief Reads an XML input file
  bool readXML(const char* fileName);
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is) = 0;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) { return false; }

  static int msgLevel; //!< Controls the amount of console output during solving

protected:
  int myPid; //!< Processor ID in parallel simulations
  int nProc; //!< Number of processors in parallel simulations

  //! \brief Certain tags needs to be parsed before others. This function takes care of this
  //! \param[in] base The base tag containing the elements to be prioritized
  //! \returns A vector with elements that was handled
  //! \details This function needs to be called in the app specific SIM classes
  //           prior parsing its app specific block. 
  std::vector<const TiXmlElement*> handlePriorityTags(const TiXmlElement* base);
};

#endif
