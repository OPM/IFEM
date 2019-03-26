// $Id$
//==============================================================================
//!
//! \file XMLInputBase.C
//!
//! \date Jul 16 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for XML input parsing functionality.
//!
//==============================================================================

#include "XMLInputBase.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <algorithm>


/*!
  \brief Helper method to load an XML file and print error message if failure.
*/

static bool loadXMLfile (TiXmlDocument& doc, const char* fileName)
{
  if (doc.LoadFile(fileName))
    return true;

  std::cerr <<" *** SIMadmin::read: Failed to load \""<< fileName
            <<"\".\n\tError at line "<< doc.ErrorRow() <<": "
            << doc.ErrorDesc() << std::endl;
  return false;
}


bool XMLInputBase::injectIncludeFiles (TiXmlElement* tag, bool verbose) const
{
  static int nLevels = 0;
  std::string spaces(2*(nLevels++),' ');
  bool foundIncludes = false, status = true;
  TiXmlElement* elem = tag->FirstChildElement();
  for (; elem && status; elem = elem->NextSiblingElement())
    if (strcasecmp(elem->Value(),"include"))
      status = this->injectIncludeFiles(elem,verbose);
    else if (elem->FirstChild() && elem->FirstChild()->Value()) {
      TiXmlDocument doc;
      if ((status = loadXMLfile(doc,elem->FirstChild()->Value()))) {
        if (verbose)
          IFEM::cout << spaces <<"Loaded included file "
                     << elem->FirstChild()->Value() << std::endl;
        TiXmlElement* e2 = doc.RootElement();
        TiXmlNode* n2 = elem = tag->ReplaceChild(elem,*e2)->ToElement();
        for (e2 = e2->NextSiblingElement(); e2; e2 = e2->NextSiblingElement())
          n2 = tag->InsertAfterChild(n2,*e2);
        foundIncludes = true;
      }
    }

  if (foundIncludes && status)
    status = this->injectIncludeFiles(tag,verbose);

  --nLevels;

  return status;
}


bool XMLInputBase::readXML (const char* fileName, bool verbose)
{
  TiXmlDocument doc;
  if (!loadXMLfile(doc,fileName))
    return false;

  const TiXmlElement* tag = doc.RootElement();
  if (!tag || strcmp(tag->Value(),"simulation")) {
    std::cerr <<" *** SIMadmin::read: Malformatted input file \""<< fileName
              <<"\"."<< std::endl;
    return false;
  }

  if (verbose)
    IFEM::cout <<"\nParsing input file "<< fileName << std::endl;

  if (!this->injectIncludeFiles(const_cast<TiXmlElement*>(tag),verbose))
    return false;

#ifdef SP_DEBUG
  if (verbose) {
    std::cout <<"\nHere is the input-file content:"<< std::endl;
    doc.Print();
  }
#endif

  std::vector<const TiXmlElement*> parsed;
  if (!this->handlePriorityTags(doc.RootElement(),parsed,verbose))
    return false;

  for (tag = tag->FirstChildElement(); tag; tag = tag->NextSiblingElement())
    if (std::find(parsed.begin(),parsed.end(),tag) == parsed.end()) {
      if (verbose)
        IFEM::cout <<"\nParsing <"<< tag->Value() <<">"<< std::endl;
      if (!this->parse(tag)) {
        std::cerr <<" *** SIMadmin::read: Failure occured while parsing \""
                  << tag->Value() <<"\"."<< std::endl;
        return false;
      }
    }

  if (verbose)
    IFEM::cout <<"\nParsing input file succeeded."<< std::endl;

  return true;
}


bool XMLInputBase::handlePriorityTags (const TiXmlElement* base,
                                       std::vector<const TiXmlElement*>& parsed,
                                       bool verbose)
{
  const char** q = this->getPrioritizedTags();
  if (!q) return true; // No prioritized tags defined

  while (*q) {
    const TiXmlElement* elem = base->FirstChildElement(*(q++));
    if (elem) {
      if (verbose)
        IFEM::cout <<"\nParsing <"<< elem->Value() <<">"<< std::endl;
      if (!this->parse(elem)) {
        std::cerr <<" *** SIMadmin::read: Failure occured while parsing \""
                  << elem->Value() <<"\"."<< std::endl;
        return false;
      }
      parsed.push_back(elem);
    }
  }

  return true;
}
