// $Id$
//==============================================================================
//!
//! \file XMLInputBase.h
//!
//! \date Jul 16 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for xml input parsing functionality.
//!
//==============================================================================

#ifndef _XML_INPUT_BASE_H
#define _XML_INPUT_BASE_H

#include <vector>

class TiXmlElement;


/*! \brief Base class for XML based input file parsing.
 *  \details This is inherited by SIMinput for input parsing handling,
 *           and is also used in applications for pre-parsing the input file.
 */
class XMLInputBase
{
public:
  //! \brief Reads an XML input file.
  //! \param fileName File to read
  //! \param verbose True to print the tags being parsed to output
  bool readXML(const char* fileName, bool verbose = true);

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) = 0;

  //! \brief Recursive helper method for processing the \a include XML-tags.
  void injectIncludeFiles(TiXmlElement* tag) const;

  //! \brief Handles the parsing order for certain XML-tags.
  //! \param[in] base The base tag containing the elements to be prioritized
  //! \param[out] parsed Vector of XML-elements that was parsed
  //! \param verbose True to print the tags being parsed to output
  //!
  //! \details Certain tags need to be parsed before others. This method takes
  //! care of this. It is called by the \a readXML method in order to read the
  //! top level tags in the required order. It can also be called by the
  //! application-specific SIM class prior to parsing its data blocks.
  //! In that case the \a getPrioritizedTags method should be reimplemented
  //! by the sub-class to take care of the application-specific tags.
  bool handlePriorityTags(const TiXmlElement* base,
                          std::vector<const TiXmlElement*>& parsed,
                          bool verbose);

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const { return nullptr; }
};

#endif
