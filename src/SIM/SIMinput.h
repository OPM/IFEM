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

class TiXmlElement;


/*!
  \brief Base class for NURBS-based FEM simulators with input file parsing.
*/

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
  //! care of this. It needs to be called in the application-specific SIM class
  //! prior to parsing its application-specific data block.
  bool handlePriorityTags(const TiXmlElement* base,
			  std::vector<const TiXmlElement*>& parsed);
  //! \brief Handles the parsing order for certain XML-tags.
  //! \param[in] base The base tag containing the elements to be prioritized
  //! \returns A vector with elements that was handled
  //!
  //! TODO: Remove this method - used the above instead.
  std::vector<const TiXmlElement*> handlePriorityTags(const TiXmlElement* base)
  {
    std::vector<const TiXmlElement*> parsed;
    handlePriorityTags(base,parsed);
    return parsed;
  }

public:
  //! \brief Parses a data section from an input stream.
  virtual bool parse(char* keyWord, std::istream& is) = 0;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) = 0;

  static int msgLevel; //!< Controls the amount of console output during solving

protected:
  int myPid; //!< Processor ID in parallel simulations
  int nProc; //!< Number of processors in parallel simulations
};

#endif
