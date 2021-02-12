// $Id$
//==============================================================================
//!
//! \file XMLInputBase.h
//!
//! \date Jul 16 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for XML input parsing functionality.
//!
//==============================================================================

#ifndef _XML_INPUT_BASE_H
#define _XML_INPUT_BASE_H

#include <vector>

class TiXmlElement;


/*!
  \brief Base class for XML based input file parsing.
  \details This class is inherited by SIMadmin for input parsing handling,
  and can also be used by applications for pre-parsing of the input file.
*/

class XMLInputBase
{
public:
  //! \brief Reads an XML input file.
  //! \param[in] fileName File to read
  //! \param[in] verbose If \e true, print the tags being parsed
  bool readXML(const char* fileName, bool verbose = true);

  //! \brief Loads data from the XML-formatted text string.
  //! \details This method is a convenience offered for unit testing only.
  bool loadXML(const char* xml);

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem) = 0;

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const { return nullptr; }

private:
  //! \brief Handles the parsing order for certain XML-tags.
  //! \param[in] base The base tag containing the elements to be prioritized
  //! \param[out] parsed Vector of XML-elements that was parsed
  //! \param[in] verbose If \e true, print the tags being parsed
  //!
  //! \details Certain tags need to be parsed before others. This method takes
  //! care of this. It is called by the readXML() method in order to read the
  //! top level tags in the required order. It can also be called by the
  //! application-specific %SIM class prior to parsing its data blocks.
  //! In that case the getPrioritizedTags() method should be reimplemented
  //! by the sub-class to take care of the application-specific tags.
  bool handlePriorityTags(const TiXmlElement* base,
                          std::vector<const TiXmlElement*>& parsed,
                          bool verbose);

  //! \brief Recursive helper method for processing the \a include XML-tags.
  bool injectIncludeFiles(TiXmlElement* tag, bool verbose) const;
};

#endif
